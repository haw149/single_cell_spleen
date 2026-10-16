import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import data_access as data


class DataAccessTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name).resolve()
        self.old_root = data.ROOT
        data.ROOT = self.root
        self.safe_patch = patch.object(data, "safe_path", lambda relative: self._safe(relative))
        self.safe_patch.start()
        self.matrix = np.array([[0, 2, 3], [0, 0, 0], [4, 0, 5], [0, 0, 0], [6, 7, 0]], dtype=float)
        obs = pd.DataFrame({"condition": pd.Categorical(["A", "B", None, "A", "B"]),
                            "nullable": pd.array([True, None, False, True, False], dtype="boolean")}, index=[f"c{i}" for i in range(5)])
        for encoding, matrix in [("csr", sparse.csr_matrix(self.matrix)), ("csc", sparse.csc_matrix(self.matrix)), ("dense", self.matrix)]:
            obj = ad.AnnData(matrix, obs=obs.copy(), var=pd.DataFrame(index=["g0", "g1", "g2"]))
            obj.obsm["X_umap"] = np.arange(10).reshape(5,2)
            obj.layers["counts"] = sparse.csr_matrix(self.matrix*2)
            obj.raw = ad.AnnData(sparse.csr_matrix(np.c_[self.matrix, [1,2,3,4,5]]), obs=obs.copy(), var=pd.DataFrame(index=["g0", "g1", "g2", "raw_only"]))
            obj.uns["test_table"] = pd.DataFrame({"gene": ["g1", "g2"], "score": [2., 3.]})
            obj.write_h5ad(self.root / f"{encoding}.h5ad")

    def _safe(self, relative):
        path = (self.root / relative).resolve()
        if not path.is_relative_to(self.root) or not path.is_file(): raise ValueError("outside root")
        return path

    def tearDown(self):
        data.clear_cache()
        self.safe_patch.stop()
        data.ROOT = self.old_root
        self.temp.cleanup()

    def test_expression_matches_dense_for_all_storage_types_and_empty_rows(self):
        for encoding in ["csr", "csc", "dense"]:
            for j in range(3):
                np.testing.assert_allclose(data.read_expression(f"{encoding}.h5ad", f"g{j}"), self.matrix[:,j])

    def test_layers_and_raw_use_correct_gene_axis(self):
        np.testing.assert_allclose(data.read_expression("csr.h5ad", "g1", "layers/counts"), self.matrix[:,1]*2)
        np.testing.assert_allclose(data.read_expression("csr.h5ad", "raw_only", "raw/X"), [1,2,3,4,5])
        with self.assertRaises(ValueError): data.read_expression("csr.h5ad", "raw_only", "X")

    def test_matrix_preview_matches_values_and_dimensions(self):
        for encoding in ["csr", "csc", "dense"]:
            preview = data.matrix_preview(f"{encoding}.h5ad", n_cells=3, n_genes=2)
            np.testing.assert_allclose(preview.iloc[:,1:].to_numpy(), self.matrix[:3,:2])
            self.assertEqual(preview.cell_id.tolist(), ["c0", "c1", "c2"])

    def test_10x_orientation_and_preview(self):
        matrix = sparse.csc_matrix(self.matrix.T)
        with h5py.File(self.root / "counts.h5", "w") as f:
            g=f.create_group("matrix")
            for key,value in dict(data=matrix.data,indices=matrix.indices,indptr=matrix.indptr,shape=matrix.shape).items(): g[key]=value
            g["barcodes"]=np.array([f"c{i}".encode() for i in range(5)])
            g["features/name"]=np.array([b"g0",b"g1",b"g2"])
        np.testing.assert_allclose(data.matrix_preview("counts.h5").iloc[:,1:].to_numpy(), self.matrix)
        self.assertEqual(len(data.load_table("counts.h5", "features")),3)

    def test_metadata_categorical_missing_values_and_embedded_tables(self):
        obs=data.load_obs("csr.h5ad")
        self.assertTrue(pd.isna(obs.condition.iloc[2]))
        self.assertTrue(pd.isna(obs.nullable.iloc[1]))
        self.assertEqual(data.describe_h5ad("csr.h5ad")["n_obs"],5)
        self.assertIn("uns/test_table",data.table_choices("csr.h5ad"))
        self.assertEqual(len(data.load_table("csr.h5ad","uns/test_table")),2)

    def test_sparse_gene_read_across_multiple_chunks(self):
        expected=np.tile(self.matrix,(250,1))
        obj=ad.AnnData(sparse.csr_matrix(expected),var=pd.DataFrame(index=["a","b","c"]))
        obj.write_h5ad(self.root/"large.h5ad")
        np.testing.assert_allclose(data.read_expression("large.h5ad","b"),expected[:,1])

    def test_scan_skips_hidden_app_and_external_symlinks(self):
        (self.root/".secret.csv").write_text("x\n1\n")
        (self.root/"shiny_app").mkdir(); (self.root/"shiny_app"/"internal.csv").write_text("x\n1\n")
        (self.root/"outside.csv").symlink_to(Path(__file__).resolve())
        paths=data.discover_files(self.root).path.tolist()
        self.assertEqual(len(paths),3)

    def test_original_safe_path_rejects_path_traversal(self):
        self.safe_patch.stop()
        try:
            with self.assertRaises(ValueError): data.safe_path("../not-here",self.root)
            self.assertEqual(data.safe_path("csr.h5ad",self.root),self.root/"csr.h5ad")
        finally: self.safe_patch.start()

    def test_table_cache_refreshes_after_save(self):
        path=self.root/"table.csv"; path.write_text("value\n1\n")
        self.assertEqual(len(data.load_table("table.csv")),1)
        path.write_text("value\n2\n3\n")
        self.assertEqual(data.load_table("table.csv").value.tolist(),[2,3])

    def test_table_limit_and_search_are_explicit(self):
        (self.root/"table.csv").write_text("gene,value\na[1],2\nb,3\nc,4\n")
        with patch.object(data,"MAX_ROWS",2):
            df=data.load_table("table.csv")
            self.assertEqual(len(df),2); self.assertIn("Preview limited",df.attrs["note"])
        self.assertEqual(data.filter_table(df,"a[1]").gene.tolist(),["a[1]"])
        self.assertEqual(data.filter_table(df,columns="value").columns.tolist(),["gene","value"])

    def test_notebook_extracts_safe_snapshot_and_preserves_cell_provenance(self):
        notebook={"cells":[{"source":["df"],"execution_count":7,"outputs":[{"data":{"text/html":["<script>alert('x')</script><table><tr><th></th><th>gene</th></tr><tr><td>0</td><td>&lt;b&gt;</td></tr></table>"],"text/plain":["truncated display"]}}]}]}
        (self.root/"notebook.ipynb").write_text(json.dumps(notebook))
        results=data.notebook_outputs("notebook.ipynb")
        self.assertEqual(results[0]["table"].gene.iloc[0],"<b>")
        self.assertEqual(results[0]["cell"],1)
        self.assertEqual(results[0]["execution"],7)
        self.assertNotIn("script",results[0]["table"].to_csv())
        self.assertEqual(results[1]["text"],"truncated display")


if __name__ == "__main__": unittest.main()
