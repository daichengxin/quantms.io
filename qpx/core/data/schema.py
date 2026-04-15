"""Schema base classes for QPX data structures.

ViewSchema instances are loaded from YAML specs at import time.
FieldDef is a simple container for a resolved Arrow field.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pyarrow as pa

# ---------------------------------------------------------------------------
# Validation result types
# ---------------------------------------------------------------------------


@dataclass
class ValidationIssue:
    """A single validation issue found during schema validation."""

    structure: str
    check: str  # "missing_column" | "type_mismatch" | "null_values" | "duplicate_pk"
    severity: str  # "error" | "warning"
    column: str | None
    message: str


@dataclass
class ValidationResult:
    """Result of validating a structure against its schema."""

    structure: str
    issues: list[ValidationIssue] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        """True if no errors (warnings are acceptable)."""
        return not any(i.severity == "error" for i in self.issues)

    @property
    def errors(self) -> list[ValidationIssue]:
        return [i for i in self.issues if i.severity == "error"]

    @property
    def warnings(self) -> list[ValidationIssue]:
        return [i for i in self.issues if i.severity == "warning"]

    @property
    def summary(self) -> str:
        n_err = len(self.errors)
        n_warn = len(self.warnings)
        status = "VALID" if self.is_valid else "INVALID"
        parts = []
        if n_err:
            parts.append(f"{n_err} error{'s' if n_err != 1 else ''}")
        if n_warn:
            parts.append(f"{n_warn} warning{'s' if n_warn != 1 else ''}")
        detail = f" ({', '.join(parts)})" if parts else ""
        return f"{self.structure}: {status}{detail}"


# ---------------------------------------------------------------------------
# Field and schema definitions
# ---------------------------------------------------------------------------


class FieldDef:
    """A resolved schema field definition."""

    __slots__ = ("name", "arrow_type", "nullable", "optional", "doc")

    def __init__(
        self,
        name: str,
        arrow_type: pa.DataType,
        *,
        nullable: bool = True,
        optional: bool = False,
        doc: str = "",
    ):
        self.name = name
        self.arrow_type = arrow_type
        self.nullable = nullable
        self.optional = optional
        self.doc = doc

    def to_arrow_field(self) -> pa.Field:
        metadata = {"doc": self.doc} if self.doc else None
        return pa.field(self.name, self.arrow_type, nullable=self.nullable, metadata=metadata)


def _types_compatible(expected: pa.DataType, actual: pa.DataType) -> bool:
    """Check type compatibility, tolerating nullability differences in nested structs.

    Parquet does not preserve inner-struct field nullability, so we compare
    types structurally while ignoring the ``nullable`` flag on struct children.
    """
    if expected == actual:
        return True

    # list<X> — compare value types recursively
    if isinstance(expected, pa.ListType) and isinstance(actual, pa.ListType):
        return _types_compatible(expected.value_type, actual.value_type)

    # map<K, V>
    if isinstance(expected, pa.MapType) and isinstance(actual, pa.MapType):
        return _types_compatible(expected.key_type, actual.key_type) and _types_compatible(expected.item_type, actual.item_type)

    # struct — compare field names, types (ignoring nullable)
    if isinstance(expected, pa.StructType) and isinstance(actual, pa.StructType):
        if expected.num_fields != actual.num_fields:
            return False
        for i in range(expected.num_fields):
            ef = expected.field(i)
            af = actual.field(i)
            if ef.name != af.name:
                return False
            if not _types_compatible(ef.type, af.type):
                return False
        return True

    return False


class ViewSchema:
    """Schema for a QPX view, loaded from YAML."""

    def __init__(
        self,
        view_name: str,
        file_type: str,
        primary_key: list[str],
        fields: dict[str, FieldDef],
        doc: str = "",
        extra_columns: bool = False,
    ):
        self._view_name = view_name
        self._file_type = file_type
        self._primary_key = tuple(primary_key)
        self._fields = fields
        self._doc = doc
        self._extra_columns = extra_columns
        self._arrow_schema_cache: pa.Schema | None = None

    @property
    def __name__(self) -> str:
        """For compatibility with parametrized tests and repr."""
        parts = self._view_name.split("_")
        return "".join(p.capitalize() for p in parts) + "Schema"

    def __repr__(self) -> str:
        return f"ViewSchema({self._view_name!r}, fields={len(self._fields)})"

    def get_arrow_schema(self) -> pa.Schema:
        """Generate PyArrow schema from field definitions."""
        if self._arrow_schema_cache is None:
            arrow_fields = [fdef.to_arrow_field() for fdef in self._fields.values()]
            self._arrow_schema_cache = pa.schema(arrow_fields)
        return self._arrow_schema_cache

    def get_extended_arrow_schema(self, extra_column_names: list[str]) -> pa.Schema:
        """Return the base schema extended with additional string columns.

        Only allowed when ``extra_columns=True`` in the YAML spec.
        Extra columns are appended as nullable ``pa.string()`` fields.
        """
        if not self._extra_columns:
            raise ValueError(f"Schema '{self._view_name}' does not allow extra columns")
        base = self.get_arrow_schema()
        known = set(base.names)
        extra_fields = [pa.field(name, pa.string(), nullable=True) for name in extra_column_names if name not in known]
        if not extra_fields:
            return base
        return pa.schema(list(base) + extra_fields)

    def _iter_fields(self):
        """Yield (name, FieldDef) pairs in definition order."""
        yield from self._fields.items()

    # --- Validation --------------------------------------------------------

    def validate(self, table: pa.Table) -> list[str]:
        """Validate an Arrow table against this schema. Returns list of error strings.

        For structured results, use validate_full() instead.
        """
        result = self.validate_full(table)
        return [issue.message for issue in result.errors]

    def validate_full(self, table: pa.Table) -> ValidationResult:
        """Full schema validation returning structured results.

        Checks:
            1. Required columns are present
            2. Column types match the schema
            3. Non-nullable columns contain no null values
            4. Primary key is unique
        """
        result = ValidationResult(structure=self._view_name)
        expected = self.get_arrow_schema()

        # Check 1 & 2: Column presence and type matching
        for expected_field in expected:
            if expected_field.name not in table.schema.names:
                if not self._is_optional(expected_field.name):
                    result.issues.append(
                        ValidationIssue(
                            structure=self._view_name,
                            check="missing_column",
                            severity="error",
                            column=expected_field.name,
                            message=f"Missing required column: {expected_field.name}",
                        )
                    )
                continue
            actual = table.schema.field(expected_field.name)
            if not _types_compatible(expected_field.type, actual.type):
                result.issues.append(
                    ValidationIssue(
                        structure=self._view_name,
                        check="type_mismatch",
                        severity="error",
                        column=expected_field.name,
                        message=(f"Column '{expected_field.name}': expected {expected_field.type}, got {actual.type}"),
                    )
                )

        # Check 3: Null values in non-nullable columns
        for field_name, fdef in self._fields.items():
            if field_name not in table.schema.names:
                continue
            if not fdef.nullable:
                col = table.column(field_name)
                null_count = col.null_count
                if null_count > 0:
                    result.issues.append(
                        ValidationIssue(
                            structure=self._view_name,
                            check="null_values",
                            severity="warning",
                            column=field_name,
                            message=(
                                f"Column '{field_name}' is non-nullable but contains "
                                f"{null_count} null value(s) out of {len(col)} rows"
                            ),
                        )
                    )

        # Check 4: Primary key uniqueness
        pk_cols = [c for c in self._primary_key if c in table.schema.names]
        if pk_cols and len(pk_cols) == len(self._primary_key) and len(table) > 0:
            # Normalize list-typed PK columns to joined strings so group_by works
            pk_arrays: dict[str, pa.Array] = {}
            for c in pk_cols:
                col = table.column(c)
                if pa.types.is_list(col.type):
                    joined = ["_".join(str(x) for x in v) if v is not None else None for v in col.to_pylist()]
                    pk_arrays[c] = pa.array(joined, type=pa.string())
                else:
                    pk_arrays[c] = col
            pk_table = pa.table(pk_arrays)
            n_total = len(pk_table)
            n_unique = pk_table.group_by(pk_cols).aggregate([]).num_rows
            if n_unique < n_total:
                n_dupes = n_total - n_unique
                result.issues.append(
                    ValidationIssue(
                        structure=self._view_name,
                        check="duplicate_pk",
                        severity="warning",
                        column=None,
                        message=(
                            f"Primary key ({', '.join(self._primary_key)}) has "
                            f"{n_dupes} duplicate row(s) "
                            f"({n_unique} unique out of {n_total})"
                        ),
                    )
                )

        return result

    def _is_optional(self, field_name: str) -> bool:
        fdef = self._fields.get(field_name)
        return fdef is not None and fdef.optional
