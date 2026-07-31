"""
Shared quantification channel-label resolution for QuantMS/OpenMS QPX.

Both QPX paths for QuantMS output — the mzTab+MSstats converter
(:mod:`qpx.converters.quantms`) and the enrichment of OpenMS ``-out_qpx``
(:mod:`qpx.converters.openms`) — must put the *same* canonical reporter labels
into ``intensities[].label`` so TMT/iTRAQ/LFQ results are consistent and join
to ``sample``/``run`` the same way DIA-NN (quantmsdiann) output does.

The canonical vocabulary is the shared sdrf-pipelines ``channel_map`` (single
source of truth for label spellings, incl. the TMT10 ch10 = ``TMT131`` vs
TMT11 ch10 = ``TMT131N`` distinction). The plex is resolved, most-solid first,
from the SDRF ``comment[label]`` set (ground truth) and only then from the
channel indices present in the data.
"""

from __future__ import annotations

import re
from typing import Iterable, Optional

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from defusedxml.common import DefusedXmlException
from defusedxml.ElementTree import ParseError, iterparse
from sdrf_pipelines.converters.channel_map import CHANNEL_LABELS
from sdrf_pipelines.converters.channel_map import normalize_label as _normalize_label
from sdrf_pipelines.converters.openms.utils import infer_itraqplex, infer_tmtplex

from qpx.writers.base import parquet_write_options

__all__ = [
    "read_sdrf_labels",
    "experiment_type_from_labels",
    "resolve_channel_labels",
    "channel_labels_from_consensusxml",
    "relabel_intensities_parquet",
    "normalize_label",
]

_FAMILY_PREFIX = {"TMT": "tmt", "ITRAQ": "itraq", "iTRAQ": "itraq"}
_CANONICAL_LABELS = {label.casefold(): label for channel_labels in CHANNEL_LABELS.values() for label in channel_labels.values()}
_ONTOLOGY_NAME = re.compile(r"(?:^|;)\s*NT\s*=\s*([^;]+)", re.IGNORECASE)


def normalize_label(label: str) -> str:
    """Return the canonical spelling of a known SDRF channel label."""
    text = label.strip()
    match = _ONTOLOGY_NAME.search(text)
    if match:
        text = match.group(1).strip()
    normalized = _normalize_label(text)
    return _CANONICAL_LABELS.get(normalized.casefold(), normalized)


def read_sdrf_labels(sdrf_path: Optional[str]) -> Optional[set[str]]:
    """
    Return the distinct ``comment[label]`` values declared in the SDRF.

    This is the experiment's ground-truth channel set. Returns ``None`` when no
    SDRF is available or it has no label column.
    """
    if not sdrf_path:
        return None
    try:
        df = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    except (OSError, pd.errors.EmptyDataError, pd.errors.ParserError):
        # Unreadable, malformed, or otherwise invalid SDRF -> no ground truth;
        # fall back to count-based inference rather than aborting the conversion.
        return None
    col = next((c for c in df.columns if c.strip().lower() == "comment[label]"), None)
    if col is None:
        return None
    labels = {str(v).strip() for v in df[col].dropna() if str(v).strip()}
    return labels or None


def experiment_type_from_labels(sdrf_labels: Optional[set[str]]) -> str:
    """Classify the labeling type (``"TMT"`` / ``"iTRAQ"`` / ``"LFQ"``) from SDRF labels."""
    if not sdrf_labels:
        return "LFQ"
    upper = {label.upper() for label in sdrf_labels}
    if any("TMT" in label for label in upper):
        return "TMT"
    if any("ITRAQ" in label for label in upper):
        return "iTRAQ"
    return "LFQ"


def resolve_channel_labels(
    experiment_type: str,
    sdrf_labels: Optional[set[str]] = None,
    channel_indices: Optional[Iterable] = None,
) -> dict[int, str]:
    """
    Resolve the ``{1-based channel index -> canonical label}`` map for a run.

    Resolution order (most-solid first):

    1. From the SDRF ``comment[label]`` set via sdrf-pipelines
       ``infer_tmtplex`` / ``infer_itraqplex`` — robust even when a channel is
       empty in the data.
    2. Fallback: the smallest plex of the reagent family whose channel count
       covers the highest index present in ``channel_indices``.

    Returns ``{}`` for non-isobaric or unresolvable input (caller falls back to
    the raw value / ``"LFQ"``).
    """
    prefix = _FAMILY_PREFIX.get(experiment_type)
    if prefix is None:
        return {}

    # 1. SDRF-declared plex (ground truth).
    if sdrf_labels:
        canonical = {normalize_label(label) for label in sdrf_labels if label}
        infer = infer_tmtplex if prefix == "tmt" else infer_itraqplex
        try:
            return CHANNEL_LABELS[infer(canonical)]
        except (ValueError, KeyError):
            pass

    # 2. Fallback: infer plex from the channel indices present in the data.
    # NB: channel_indices may be a numpy array (df[col].values), so test for
    # None explicitly — `channel_indices or ()` raises on array truthiness.
    indices: set[int] = set()
    for value in () if channel_indices is None else channel_indices:
        try:
            indices.add(int(float(value)))
        except (TypeError, ValueError):
            continue
    if not indices:
        return {}
    max_index = max(indices)

    family = {plex: labels for plex, labels in CHANNEL_LABELS.items() if plex.startswith(prefix)}
    for plex in sorted(family, key=lambda p: len(family[p])):
        if len(family[plex]) >= max_index:
            return family[plex]
    return {}


def channel_labels_from_consensusxml(
    consensusxml_path: str,
    experiment_type: str,
    sdrf_labels: Optional[set[str]] = None,
) -> dict[int, str]:
    """
    Resolve ``{1-based channel index -> canonical label}`` from a consensusXML.

    The OpenMS ConsensusMap ``ColumnHeaders`` give the **authoritative** channel
    count and order (map-index), independent of which channels happen to carry
    data — so this is more robust than inferring the plex from the data alone.
    The canonical spelling still comes from the SDRF-declared plex + channel_map;
    a ColumnHeader ``label`` is used only as a fallback for a position the plex
    map does not cover. Only the leading ``mapList`` is parsed, so even very
    large consensusXML files are handled with bounded memory. Returns ``{}``
    when the file is unavailable or malformed.
    """
    headers: dict[int, str] = {}
    try:
        in_map_list = False
        for event, element in iterparse(consensusxml_path, events=("start", "end")):
            tag = element.tag.rsplit("}", 1)[-1]
            if event == "start" and tag == "mapList":
                in_map_list = True
            elif event == "start" and in_map_list and tag == "map":
                map_index = int(element.attrib["id"])
                headers[map_index] = element.attrib.get("label", "").strip()
            elif event == "end" and tag == "mapList":
                break
            if event == "end":
                element.clear()
    except (OSError, ParseError, DefusedXmlException, KeyError, ValueError):
        return {}
    if not headers:
        return {}

    n_channels = len(headers)
    # Plex resolved with the consensusXML channel count as the index range.
    plex = resolve_channel_labels(experiment_type, sdrf_labels, range(1, n_channels + 1))
    labels: dict[int, str] = {}
    for map_index in sorted(headers):
        position = map_index + 1  # ColumnHeaders are 0-based; qpx channels are 1-based
        raw = headers[map_index]
        labels[position] = plex.get(position) or (normalize_label(raw) if raw else str(position))
    return labels


def _relabel_entries(rows, channel_labels: dict[int, str], is_lfq: bool):
    """
    Relabel the ``label`` field of a list<struct> intensities column.

    OpenMS ``-out_qpx`` writes the *filename* into feature ``intensities[].label``
    and a bare 1-based index (``"1"``, ``"2"``, ...) into pg ``intensities[].label``.
    Recover the channel index from a numeric label, else from the array position
    (channels are emitted in index order), then map to the canonical label.
    LFQ collapses to the single ``"LFQ"`` label.
    """
    out = []
    for row in rows:
        if row is None:
            out.append(None)
            continue
        new_row = []
        for position, entry in enumerate(row):
            label = entry.get("label")
            if is_lfq:
                new_label = "LFQ"
            else:
                text = str(label) if label is not None else ""
                index = int(text) if text.isdigit() else position + 1
                new_label = channel_labels.get(index, normalize_label(text) if text else text)
            updated = dict(entry)
            updated["label"] = new_label
            new_row.append(updated)
        out.append(new_row)
    return out


def _resolve_parquet_compression(
    parquet: pq.ParquetFile,
    source_metadata: dict[bytes, bytes],
    compression: str | None,
) -> str:
    """Resolve the output codec from an override, QPX metadata, or the file."""
    if compression is not None:
        return compression

    declared = source_metadata.get(b"compression_format")
    if declared:
        return declared.decode()

    if parquet.num_row_groups:
        row_group = parquet.metadata.row_group(0)
        codecs = {row_group.column(index).compression.lower() for index in range(row_group.num_columns)}
        if len(codecs) == 1:
            codec = codecs.pop()
            return "none" if codec == "uncompressed" else codec
    return "snappy"


def relabel_intensities_parquet(
    src_path: str,
    dst_path: str,
    channel_labels: dict[int, str],
    is_lfq: bool,
    columns: tuple[str, ...] = ("intensities", "additional_intensities"),
    compression: str | None = None,
) -> None:
    """
    Rewrite ``src_path`` to ``dst_path`` with canonical channel labels.

    Streams by row group so memory stays bounded on large feature tables.
    Non-intensity columns pass through untouched. QPX compression and selective
    BYTE_STREAM_SPLIT encoding are retained; *compression* overrides the source
    codec when supplied.
    """
    parquet = pq.ParquetFile(src_path)
    source_metadata = dict(parquet.schema_arrow.metadata or {})
    output_compression = _resolve_parquet_compression(parquet, source_metadata, compression)

    if compression is not None or b"compression_format" in source_metadata:
        source_metadata[b"compression_format"] = output_compression.encode()
    output_schema = parquet.schema_arrow.with_metadata(source_metadata or None)
    writer = pq.ParquetWriter(
        dst_path,
        output_schema,
        **parquet_write_options(output_schema, output_compression),
    )
    try:
        for group in range(parquet.num_row_groups):
            table = parquet.read_row_group(group)
            for column in columns:
                if column not in table.column_names:
                    continue
                field_index = table.schema.get_field_index(column)
                original = table.column(column)
                relabeled = pa.array(
                    _relabel_entries(original.to_pylist(), channel_labels, is_lfq),
                    type=original.type,
                )
                table = table.set_column(field_index, column, relabeled)
            table = table.replace_schema_metadata(output_schema.metadata)
            writer.write_table(table)
    finally:
        writer.close()
