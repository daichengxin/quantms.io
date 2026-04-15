"""DatasetCollection tests."""

import shutil


class TestVirtualCollection:
    def test_sql_across_datasets(self, dataset_dir, tmp_path):
        import qpx

        ds2_dir = tmp_path / "ds2"
        shutil.copytree(dataset_dir, ds2_dir)

        ds1 = qpx.open_dataset(str(dataset_dir))
        ds2 = qpx.open_dataset(str(ds2_dir))

        coll = qpx.DatasetCollection([ds1, ds2])
        result = coll.sql("SELECT COUNT(*) AS cnt FROM feature_0")
        assert result.fetchone()[0] > 0
        coll.close()
        ds1.close()
        ds2.close()

    def test_sql_union_across_datasets(self, dataset_dir, tmp_path):
        import qpx

        ds2_dir = tmp_path / "ds2"
        shutil.copytree(dataset_dir, ds2_dir)

        ds1 = qpx.open_dataset(str(dataset_dir))
        ds2 = qpx.open_dataset(str(ds2_dir))

        coll = qpx.DatasetCollection([ds1, ds2])
        result = coll.sql("SELECT COUNT(*) AS cnt FROM feature_0 UNION ALL SELECT COUNT(*) FROM feature_1")
        rows = result.fetchall()
        assert len(rows) == 2
        assert rows[0][0] == rows[1][0]  # Same data, same count
        coll.close()
        ds1.close()
        ds2.close()

    def test_structure_names(self, dataset_dir, tmp_path):
        import qpx

        ds2_dir = tmp_path / "ds2"
        shutil.copytree(dataset_dir, ds2_dir)

        ds1 = qpx.open_dataset(str(dataset_dir))
        ds2 = qpx.open_dataset(str(ds2_dir))

        coll = qpx.DatasetCollection([ds1, ds2])
        names = coll.structure_names
        assert 0 in names
        assert 1 in names
        assert "feature" in names[0]
        assert "feature" in names[1]
        coll.close()
        ds1.close()
        ds2.close()

    def test_context_manager(self, dataset_dir, tmp_path):
        import qpx

        ds2_dir = tmp_path / "ds2"
        shutil.copytree(dataset_dir, ds2_dir)

        ds1 = qpx.open_dataset(str(dataset_dir))
        ds2 = qpx.open_dataset(str(ds2_dir))

        with qpx.DatasetCollection([ds1, ds2]) as coll:
            result = coll.sql("SELECT COUNT(*) FROM feature_0")
            assert result.fetchone()[0] > 0

        ds1.close()
        ds2.close()


class TestPhysicalMerge:
    def test_merge_creates_output(self, dataset_dir, tmp_path):
        import qpx

        ds2_dir = tmp_path / "ds2"
        shutil.copytree(dataset_dir, ds2_dir)
        merge_dir = tmp_path / "merged"

        ds1 = qpx.open_dataset(str(dataset_dir))
        ds2 = qpx.open_dataset(str(ds2_dir))

        coll = qpx.DatasetCollection([ds1, ds2])
        coll.merge(merge_dir, structures=["feature"])

        # Verify merged file exists and has double rows
        merged_ds = qpx.open_dataset(str(merge_dir))
        assert merged_ds.feature is not None
        assert merged_ds.feature.count() == ds1.feature.count() + ds2.feature.count()

        coll.close()
        ds1.close()
        ds2.close()
        merged_ds.close()

    def test_merge_adds_source_column(self, dataset_dir, tmp_path):
        import qpx

        ds2_dir = tmp_path / "ds2"
        shutil.copytree(dataset_dir, ds2_dir)
        merge_dir = tmp_path / "merged"

        ds1 = qpx.open_dataset(str(dataset_dir))
        ds2 = qpx.open_dataset(str(ds2_dir))

        coll = qpx.DatasetCollection([ds1, ds2])
        coll.merge(merge_dir, structures=["feature"])

        merged_ds = qpx.open_dataset(str(merge_dir))
        df = merged_ds.feature.to_df()
        assert "source_dataset" in df.columns
        assert df["source_dataset"].nunique() == 2

        coll.close()
        ds1.close()
        ds2.close()
        merged_ds.close()

    def test_merge_common_structures(self, dataset_dir, tmp_path):
        import qpx

        ds2_dir = tmp_path / "ds2"
        shutil.copytree(dataset_dir, ds2_dir)
        merge_dir = tmp_path / "merged"

        ds1 = qpx.open_dataset(str(dataset_dir))
        ds2 = qpx.open_dataset(str(ds2_dir))

        coll = qpx.DatasetCollection([ds1, ds2])
        # Merge all common structures
        coll.merge(merge_dir)

        merged_ds = qpx.open_dataset(str(merge_dir))
        # At minimum, feature should be present
        assert merged_ds.feature is not None

        coll.close()
        ds1.close()
        ds2.close()
        merged_ds.close()
