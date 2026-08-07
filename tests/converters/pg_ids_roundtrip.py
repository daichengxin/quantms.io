"""Shared read-back + assertions for the feature.pg_ids round-trip checks.

Used by the openms-consensus, DIA-NN and MaxQuant converter tests to verify the
feature->pg cross-reference (bigbio/qpx#266): every populated ``feature.pg_ids``
value is a real ``pg.pg_id`` whose ``pg_accessions`` equals the feature's group
and whose ``grouped_runs`` contains the feature's ``run_file_name``.
"""

from __future__ import annotations

import duckdb


def read_feature_pg(feature_path, pg_path):
    """Read back (feature rows, pg rows) for the feature.pg_ids checks.

    feature rows: ``(sequence, run_file_name, [pg_accessions...], [pg_ids...])``;
    pg rows: ``(pg_id, [pg_accessions...], [grouped_runs...], label)``.
    ``feature.pg_accessions`` is a list of structs, so its accessions are
    projected out to a plain string list to compare against the pg view's
    ``list<string>`` membership.
    """
    con = duckdb.connect()
    feat = con.execute(
        "SELECT sequence, run_file_name, list_transform(pg_accessions, x -> x.accession), pg_ids FROM read_parquet($1)",
        [str(feature_path)],
    ).fetchall()
    pg = con.execute(
        "SELECT pg_id, pg_accessions, grouped_runs, label FROM read_parquet($1)",
        [str(pg_path)],
    ).fetchall()
    con.close()
    return feat, pg


def assert_pg_ids_join_valid(feat, pg):
    """Every populated feature.pg_ids id references a real pg row whose
    pg_accessions equals the feature's (as a set — both are membership sets) and
    whose grouped_runs contains the feature's run.
    """
    pg_by_id = {pg_id: (list(accs), list(grouped_runs)) for pg_id, accs, grouped_runs, _ in pg}
    for _seq, run, accs, pg_ids in feat:
        if not pg_ids:
            continue
        feat_accs = set(accs) if accs is not None else set()
        for pg_id in pg_ids:
            assert pg_id in pg_by_id, f"feature.pg_ids references unknown pg_id {pg_id}"
            row_accs, row_runs = pg_by_id[pg_id]
            assert set(row_accs) == feat_accs, f"pg row {pg_id} accessions {row_accs} != feature {sorted(feat_accs)}"
            assert run in row_runs, f"feature run {run!r} not in pg row {pg_id} grouped_runs {row_runs}"
