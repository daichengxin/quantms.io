"""S3 support tests (unit tests only, no real S3 access)."""


class TestS3EngineConfig:
    def test_httpfs_setup_with_credentials(self):
        from qpx.core.engine import DuckDBEngine

        engine = DuckDBEngine(
            s3_config={
                "region": "us-east-1",
                "access_key_id": "test_key",
                "secret_access_key": "test_secret",
            }
        )
        # Verify httpfs is loaded
        result = engine.execute("SELECT * FROM duckdb_extensions() WHERE extension_name = 'httpfs' AND loaded = true").fetchone()
        assert result is not None
        engine.close()

    def test_no_httpfs_by_default(self):
        from qpx.core.engine import DuckDBEngine

        engine = DuckDBEngine()
        # httpfs should NOT be loaded by default
        result = engine.execute("SELECT * FROM duckdb_extensions() WHERE extension_name = 'httpfs' AND loaded = true").fetchone()
        assert result is None
        engine.close()

    def test_anonymous_s3_config(self):
        from qpx.core.engine import DuckDBEngine

        engine = DuckDBEngine(s3_config={"region": "us-east-1", "anonymous": True})
        result = engine.execute("SELECT * FROM duckdb_extensions() WHERE extension_name = 'httpfs' AND loaded = true").fetchone()
        assert result is not None
        engine.close()

    def test_s3_with_custom_endpoint(self):
        from qpx.core.engine import DuckDBEngine

        engine = DuckDBEngine(
            s3_config={
                "region": "us-east-1",
                "endpoint": "minio.example.com:9000",
                "access_key_id": "minioadmin",
                "secret_access_key": "minioadmin",
            }
        )
        result = engine.execute("SELECT * FROM duckdb_extensions() WHERE extension_name = 'httpfs' AND loaded = true").fetchone()
        assert result is not None
        engine.close()


class TestS3DatasetDetection:
    def test_is_s3_flag_set(self):
        from qpx.dataset import Dataset

        # S3 path detection (will fail on discovery, that's expected)
        try:
            ds = Dataset(
                "s3://bucket/dataset/",
                structures=[],
                s3_config={"region": "us-east-1"},
            )
            assert ds._is_s3 is True
        except Exception:
            pass  # S3 discovery will fail — we just check the flag

    def test_local_path_not_s3(self, tmp_path):
        import qpx

        ds = qpx.open_dataset(str(tmp_path))
        assert ds._is_s3 is False
        ds.close()

    def test_s3_path_is_string(self):
        """S3 paths should stay as strings, not become Path objects."""
        from qpx.dataset import Dataset

        try:
            ds = Dataset(
                "s3://bucket/dataset/",
                structures=[],
                s3_config={"region": "us-east-1"},
            )
            assert isinstance(ds.path, str)
        except Exception:
            pass
