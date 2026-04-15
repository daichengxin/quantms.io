"""Tests for qpx.core.sql — SQL safety utilities."""

import pytest

from qpx.core.sql import escape_path, sql_build, validate_identifier, validate_table

# ---------------------------------------------------------------------------
# validate_identifier
# ---------------------------------------------------------------------------


class TestValidateIdentifier:
    """validate_identifier: double-quotes safe column names."""

    def test_simple_name(self):
        assert validate_identifier("run_file_name") == '"run_file_name"'

    def test_bracket_name(self):
        assert validate_identifier("best_search_engine_score[1]") == '"best_search_engine_score[1]"'

    def test_name_with_space(self):
        assert validate_identifier("Ion Mobility") == '"Ion Mobility"'

    def test_name_with_dot(self):
        assert validate_identifier("table.column") == '"table.column"'

    def test_rejects_semicolon(self):
        with pytest.raises(ValueError, match="Unsafe SQL identifier"):
            validate_identifier("name; DROP TABLE")

    def test_rejects_single_quote(self):
        with pytest.raises(ValueError, match="Unsafe SQL identifier"):
            validate_identifier("it's")

    def test_rejects_double_dash(self):
        with pytest.raises(ValueError, match="Unsafe SQL identifier"):
            validate_identifier("col--comment")

    def test_rejects_empty(self):
        with pytest.raises(ValueError, match="Unsafe SQL identifier"):
            validate_identifier("")


# ---------------------------------------------------------------------------
# validate_table
# ---------------------------------------------------------------------------


class TestValidateTable:
    """validate_table: allows bare alphanumeric/underscore names."""

    def test_simple_table(self):
        assert validate_table("report") == "report"

    def test_underscore_prefix(self):
        assert validate_table("_psm_lookup") == "_psm_lookup"

    def test_rejects_dot(self):
        with pytest.raises(ValueError, match="Unsafe table"):
            validate_table("schema.table")

    def test_rejects_space(self):
        with pytest.raises(ValueError, match="Unsafe table"):
            validate_table("my table")

    def test_rejects_sql_injection(self):
        with pytest.raises(ValueError, match="Unsafe table"):
            validate_table("t; DROP TABLE users")

    def test_rejects_empty(self):
        with pytest.raises(ValueError, match="Unsafe table"):
            validate_table("")


# ---------------------------------------------------------------------------
# escape_path
# ---------------------------------------------------------------------------


class TestEscapePath:
    """escape_path: single-quote escaping for SQL string literals."""

    def test_no_quotes(self):
        assert escape_path("/data/file.parquet") == "/data/file.parquet"

    def test_single_quote(self):
        assert escape_path("C:/data/file's.parquet") == "C:/data/file''s.parquet"

    def test_multiple_quotes(self):
        assert escape_path("a'b'c") == "a''b''c"

    def test_backslash_path(self):
        assert escape_path("D:\\data\\file.parquet") == "D:\\data\\file.parquet"


# ---------------------------------------------------------------------------
# sql_build
# ---------------------------------------------------------------------------


class TestSqlBuild:
    """sql_build: Template-based safe SQL construction."""

    def test_simple_substitution(self):
        result = sql_build("SELECT * FROM $table", table="report")
        assert result == "SELECT * FROM report"

    def test_multiple_substitutions(self):
        result = sql_build(
            "SELECT $col FROM $table WHERE $col IS NOT NULL",
            col='"run_file_name"',
            table="msstats",
        )
        assert result == 'SELECT "run_file_name" FROM msstats WHERE "run_file_name" IS NOT NULL'

    def test_preserves_duckdb_positional(self):
        result = sql_build("SELECT * FROM $table WHERE id = $1", table="users")
        assert result == "SELECT * FROM users WHERE id = $1"

    def test_preserves_question_mark(self):
        result = sql_build("SELECT * FROM $table WHERE name = ?", table="users")
        assert result == "SELECT * FROM users WHERE name = ?"

    def test_unreplaced_var_preserved(self):
        result = sql_build("SELECT $a, $b", a="col1")
        assert result == "SELECT col1, $b"

    def test_combined_with_validators(self):
        col = validate_identifier("intensity")
        tbl = validate_table("features")
        path = escape_path("C:/my file's.parquet")
        result = sql_build(
            "CREATE VIEW $t AS SELECT $c FROM read_parquet('$p')",
            t=tbl,
            c=col,
            p=path,
        )
        assert result == """CREATE VIEW features AS SELECT "intensity" FROM read_parquet('C:/my file''s.parquet')"""
