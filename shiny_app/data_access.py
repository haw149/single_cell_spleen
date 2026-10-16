"""Read-only access to project artifacts, with bounded sparse-matrix reads."""
from __future__ import annotations

import json
import os
from datetime import datetime
from functools import lru_cache
from html.parser import HTMLParser
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from scipy import sparse

ROOT = Path(os.environ.get("EDNA_PROJECT_ROOT", Path(__file__).resolve().parents[1])).resolve()
MAX_ROWS = 100_000
SKIP_DIRS = {"shiny_app", "node_modules", "__pycache__", "renv", "venv"}


def safe_path(relative: str, root: Path = ROOT) -> Path:
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()) or not path.is_file():
        raise ValueError("Select an existing file inside the project folder.")
    return path


def discover_files(root: Path = ROOT) -> pd.DataFrame:
    rows = []
    for directory, folders, files in os.walk(root, followlinks=False):
        folders[:] = sorted(d for d in folders if not d.startswith(".") and d not in SKIP_DIRS)
        for name in sorted(files):
            if name.startswith("."):
                continue
            path = Path(directory) / name
            if not path.resolve().is_relative_to(root.resolve()):
                continue
            suffix = path.suffix.lower()
            kind = ({".h5ad": "AnnData", ".h5": "HDF5", ".hdf5": "HDF5",
                     ".ipynb": "Notebook", ".csv": "Table", ".tsv": "Table",
                     ".parquet": "Table", ".xlsx": "Table", ".txt": "Text",
                     ".json": "Text", ".png": "Figure", ".jpg": "Figure",
                     ".jpeg": "Figure", ".svg": "Figure", ".pdf": "Figure"}.get(suffix, "Other"))
            try:
                stat = path.stat()
            except OSError:
                continue  # A producer may rename a file between directory scan and stat.
            rows.append(dict(path=str(path.relative_to(root)), kind=kind, size_bytes=stat.st_size,
                             modified=datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M"),
                             mtime_ns=stat.st_mtime_ns))
    return pd.DataFrame(rows, columns=["path", "kind", "size_bytes", "modified", "mtime_ns"])


def stamp(relative: str) -> tuple:
    path = safe_path(relative)
    stat = path.stat()
    return str(path), stat.st_mtime_ns, stat.st_size


def _strings(values):
    return np.array([x.decode("utf-8", errors="replace") if isinstance(x, bytes) else x for x in values])


def _shape(node):
    return tuple(node.shape or ()) if isinstance(node, h5py.Dataset) else tuple(node.attrs.get("shape", ()))


def _frame(group, index_name):
    df = read_elem(group)
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)
    df = df.copy()
    if index_name in df.columns:
        index_name = "index_" + index_name
    df.insert(0, index_name, df.index.astype(str))
    return df.reset_index(drop=True)


@lru_cache(maxsize=12)
def _describe(token):
    with h5py.File(token[0], "r") as f:
        shape = _shape(f["X"])
        layers = ["X"] + ["layers/" + key for key in f.get("layers", {})]
        if "raw/X" in f:
            layers.append("raw/X")
        tables = {"obs": "Cell metadata · obs", "var": "Gene metadata · var"}
        if "raw/var" in f:
            tables["raw/var"] = "Raw gene metadata · raw/var"
        matrices = []

        def inspect(name, node):
            encoding = node.attrs.get("encoding-type", "")
            shape_ = _shape(node)
            if encoding == "dataframe" and name not in tables:
                tables[name] = name
            elif isinstance(node, h5py.Dataset) and name.startswith("uns/") and node.size is not None and node.size <= MAX_ROWS:
                if node.ndim <= 2:
                    tables[name] = name
            if encoding in ("csr_matrix", "csc_matrix") or (isinstance(node, h5py.Dataset) and len(shape_) == 2):
                matrices.append({"location": name, "shape": " × ".join(map(str, shape_)),
                                 "storage": encoding or "dense", "dtype": str(node["data"].dtype if isinstance(node, h5py.Group) else node.dtype)})

        f.visititems(inspect)
        embeddings = [key for key in f.get("obsm", {}) if len(_shape(f["obsm"][key])) == 2 and _shape(f["obsm"][key])[1] >= 2]
        for key in embeddings:
            tables["obsm/" + key] = "Embedding · " + key
        return dict(n_obs=shape[0], n_vars=shape[1], embeddings=embeddings, layers=layers,
                    obs_columns=list(_strings(f["obs"].attrs["column-order"])),
                    var_columns=list(_strings(f["var"].attrs["column-order"])),
                    tables=tables, matrices=pd.DataFrame(matrices))


def describe_h5ad(relative):
    return _describe(stamp(relative))


@lru_cache(maxsize=12)
def _metadata(token, key):
    with h5py.File(token[0], "r") as f:
        return _frame(f[key], "cell_id" if key == "obs" else "gene")


def load_obs(relative):
    return _metadata(stamp(relative), "obs").copy()


def load_var(relative, raw=False):
    return _metadata(stamp(relative), "raw/var" if raw else "var").copy()


@lru_cache(maxsize=8)
def _embedding(token, key):
    with h5py.File(token[0], "r") as f:
        return np.asarray(f["obsm"][key][:, :2])


def read_embedding(relative, key):
    return _embedding(stamp(relative), key)


def gene_names(relative, layer="X"):
    return load_var(relative, layer == "raw/X")["gene"].tolist()


def _column(node, column):
    """Read one matrix column without materializing the entire matrix."""
    shape = _shape(node)
    if isinstance(node, h5py.Dataset):
        return np.asarray(node[:, column]).ravel()
    encoding = node.attrs.get("encoding-type")
    values = np.zeros(shape[0], dtype=np.float64)
    if encoding == "csc_matrix":
        lo, hi = node["indptr"][column:column + 2]
        np.add.at(values, node["indices"][lo:hi], node["data"][lo:hi])
    elif encoding == "csr_matrix":
        ptr = node["indptr"][:]
        for start in range(0, shape[0], 512):
            stop = min(start + 512, shape[0])
            lo, hi = int(ptr[start]), int(ptr[stop])
            indices = node["indices"][lo:hi]
            hits = np.flatnonzero(indices == column)
            if len(hits):
                data = node["data"][lo:hi]
                rows = np.searchsorted(ptr[start:stop + 1] - lo, hits, side="right") - 1 + start
                np.add.at(values, rows, data[hits])
    else:
        raise ValueError(f"Unsupported matrix encoding: {encoding}")
    return values


@lru_cache(maxsize=12)
def _expression(token, gene, layer):
    with h5py.File(token[0], "r") as f:
        var = f["raw/var" if layer == "raw/X" else "var"]
        names = _strings(var[var.attrs["_index"]][:]).tolist()
        if gene not in names:
            raise ValueError(f"Gene {gene!r} is absent from this layer.")
        return _column(f[layer], names.index(gene))


def read_expression(relative, gene, layer="X"):
    return _expression(stamp(relative), gene, layer)


def h5_inventory(relative):
    rows = []
    with h5py.File(safe_path(relative), "r") as f:
        def visit(name, node):
            if isinstance(node, h5py.Dataset):
                rows.append({"location": name, "shape": " × ".join(map(str, node.shape or ())) or "scalar", "dtype": str(node.dtype)})
        f.visititems(visit)
    return pd.DataFrame(rows)


def table_choices(relative):
    suffix = safe_path(relative).suffix.lower()
    if suffix == ".h5ad":
        return {**describe_h5ad(relative)["tables"], "__matrix__": "Expression matrix · first 50 cells × 20 genes", "__schema__": "Matrix inventory"}
    if suffix in (".h5", ".hdf5"):
        with h5py.File(safe_path(relative), "r") as f:
            if "matrix/features" in f and "matrix/barcodes" in f:
                return {"features": "Feature metadata", "barcodes": "Cell barcodes · first 100,000", "__matrix__": "Counts · first 50 cells × 20 genes", "__schema__": "HDF5 inventory"}
        return {"__schema__": "HDF5 inventory"}
    if suffix == ".xlsx":
        with pd.ExcelFile(safe_path(relative)) as book:
            return {name: name for name in book.sheet_names}
    return {"table": "Data table"}


def matrix_preview(relative, layer="X", n_cells=50, n_genes=20):
    with h5py.File(safe_path(relative), "r") as f:
        if "matrix" in f and "obs" not in f:
            g = f["matrix"]
            nc = min(n_cells, len(g["barcodes"]))
            ng = min(n_genes, int(g["shape"][0]))
            ptr = g["indptr"][:nc + 1]
            end = int(ptr[-1])
            block = sparse.csc_matrix((g["data"][:end], g["indices"][:end], ptr), shape=(int(g["shape"][0]), nc))[:ng].toarray().T
            cells = _strings(g["barcodes"][:nc])
            genes = _strings(g["features/name"][:ng])
        else:
            node = f[layer]
            nc, ng = min(n_cells, _shape(node)[0]), min(n_genes, _shape(node)[1])
            if isinstance(node, h5py.Dataset):
                block = node[:nc, :ng]
            elif node.attrs["encoding-type"] == "csr_matrix":
                ptr = node["indptr"][:nc + 1]
                end = int(ptr[-1])
                block = sparse.csr_matrix((node["data"][:end], node["indices"][:end], ptr), shape=(nc, _shape(node)[1]))[:, :ng].toarray()
            else:
                ptr = node["indptr"][:ng + 1]
                end = int(ptr[-1])
                block = sparse.csc_matrix((node["data"][:end], node["indices"][:end], ptr), shape=(_shape(node)[0], ng))[:nc].toarray()
            obs, var = f["obs"], f["raw/var" if layer == "raw/X" else "var"]
            cells = _strings(obs[obs.attrs["_index"]][:nc])
            genes = _strings(var[var.attrs["_index"]][:ng])
        # Keep duplicated feature symbols distinguishable in the grid.
        labels = [str(g) if list(genes).count(g) == 1 else f"{g} [{i}]" for i, g in enumerate(genes)]
        df = pd.DataFrame(block, columns=labels)
        df.insert(0, "cell_id", cells)
        df.attrs["note"] = f"Matrix preview only: first {nc:,} cells and {ng:,} genes. Values are read from {layer}."
        return df


@lru_cache(maxsize=6)
def _table(token, key):
    path = Path(token[0])
    suffix = path.suffix.lower()
    note = ""
    if suffix in (".csv", ".tsv"):
        df = pd.read_csv(path, sep="\t" if suffix == ".tsv" else ",", nrows=MAX_ROWS + 1, low_memory=False)
    elif suffix == ".xlsx":
        df = pd.read_excel(path, sheet_name=key, nrows=MAX_ROWS + 1)
    elif suffix == ".parquet":
        import pyarrow.parquet as pq
        parquet = pq.ParquetFile(path)
        batches, count = [], 0
        for batch in parquet.iter_batches(batch_size=10_000):
            batches.append(batch.to_pandas())
            count += len(batch)
            if count > MAX_ROWS:
                break
        df = pd.concat(batches, ignore_index=True) if batches else pd.DataFrame(columns=parquet.schema.names)
    elif suffix in (".h5", ".hdf5", ".h5ad"):
        with h5py.File(path, "r") as f:
            if key == "features":
                group = f["matrix/features"]
                n = int(f["matrix/shape"][0])
                df = pd.DataFrame({k: _strings(v[:MAX_ROWS + 1]) for k, v in group.items() if isinstance(v, h5py.Dataset) and v.shape == (n,)})
            elif key == "barcodes":
                df = pd.DataFrame({"cell_id": _strings(f["matrix/barcodes"][:MAX_ROWS + 1])})
            else:
                node = f[key]
                if node.attrs.get("encoding-type") == "dataframe":
                    df = _frame(node, "cell_id" if key == "obs" else "gene" if key.endswith("var") else "index")
                elif isinstance(node, h5py.Dataset):
                    value = node[()] if node.ndim == 0 else node[:MAX_ROWS + 1]
                    if node.dtype.names:
                        df = pd.DataFrame.from_records(np.atleast_1d(value))
                    elif node.ndim == 0:
                        df = pd.DataFrame({"value": [value.decode() if isinstance(value, bytes) else value]})
                    elif node.ndim == 1:
                        df = pd.DataFrame({"value": _strings(value)})
                    else:
                        df = pd.DataFrame(value, columns=[f"component_{i+1}" for i in range(value.shape[1])])
                        if key.startswith("obsm/"):
                            obs = f["obs"]
                            df.insert(0, "cell_id", _strings(obs[obs.attrs["_index"]][:len(df)]))
                    for col in df:
                        if df[col].dtype == object:
                            df[col] = df[col].map(lambda v: v.decode(errors="replace") if isinstance(v, bytes) else v)
                else:
                    raise ValueError("This object is not a tabular dataset.")
    else:
        raise ValueError("This file is available through its original download; it has no table reader.")
    if len(df) > MAX_ROWS:
        df = df.iloc[:MAX_ROWS].copy()
        note = f"Preview limited to the first {MAX_ROWS:,} rows; download the original file for all data."
    df.columns = [str(c) for c in df.columns]
    df.attrs["note"] = note
    return df


def load_table(relative, key="table"):
    if key == "__matrix__":
        return matrix_preview(relative)
    if key == "__schema__":
        return describe_h5ad(relative)["matrices"] if relative.endswith(".h5ad") else h5_inventory(relative)
    return _table(stamp(relative), key).copy()


class _TableParser(HTMLParser):
    """Extract displayed cell text only; notebook HTML is never executed."""
    def __init__(self):
        super().__init__()
        self.tables, self.rows, self.row, self.cell = [], None, None, None

    def handle_starttag(self, tag, attrs):
        if tag == "table": self.rows = []
        elif tag == "tr" and self.rows is not None: self.row = []
        elif tag in ("th", "td") and self.row is not None: self.cell = []

    def handle_data(self, data):
        if self.cell is not None: self.cell.append(data)

    def handle_endtag(self, tag):
        if tag in ("th", "td") and self.cell is not None:
            self.row.append("".join(self.cell).strip()); self.cell = None
        elif tag == "tr" and self.row is not None:
            self.rows.append(self.row); self.row = None
        elif tag == "table" and self.rows is not None:
            self.tables.append(self.rows); self.rows = None


def _join(value):
    return "".join(value) if isinstance(value, list) else str(value or "")


@lru_cache(maxsize=4)
def _notebook(token):
    notebook = json.loads(Path(token[0]).read_text())
    results = []
    for i, cell in enumerate(notebook.get("cells", [])):
        source = _join(cell.get("source", ""))
        for j, output in enumerate(cell.get("outputs", [])):
            data = output.get("data", {})
            entry = dict(cell=i + 1, output=j + 1, execution=cell.get("execution_count"), source=source)
            if "text/html" in data:
                parser = _TableParser(); parser.feed(_join(data["text/html"]))
                for k, rows in enumerate(parser.tables):
                    if not rows: continue
                    width = max(map(len, rows))
                    header = rows[0] + [""] * (width - len(rows[0]))
                    header = [h or f"index_{c+1}" for c, h in enumerate(header)]
                    header = [h if header.count(h) == 1 else f"{h}_{c+1}" for c,h in enumerate(header)]
                    df = pd.DataFrame([r + [""] * (width-len(r)) for r in rows[1:]], columns=header)
                    results.append({**entry, "kind": "Table snapshot", "table": df, "label": f"Cell {i+1} · table {k+1}"})
            if "image/png" in data:
                results.append({**entry, "kind": "Figure snapshot", "png": _join(data["image/png"]), "label": f"Cell {i+1} · figure {j+1}"})
            text = _join(output.get("text", data.get("text/plain", "")))
            if output.get("output_type") == "error":
                text = f"{output.get('ename', 'Error')}: {output.get('evalue', '')}"
            if text.strip():
                results.append({**entry, "kind": "Text snapshot", "text": text, "label": f"Cell {i+1} · text {j+1}"})
    return results


def notebook_outputs(relative):
    return _notebook(stamp(relative))


def filter_table(df, query="", columns="", limit=60):
    """Search loaded rows before limiting visible columns, using literal text."""
    if query.strip():
        query = query.strip()
        mask = pd.Series(False, index=df.index)
        # Wide reference matrices are searched row-wise to avoid 30,000 Python loops.
        if len(df.columns) > 200:
            mask = df.apply(lambda row: row.astype(str).str.contains(query, case=False, regex=False, na=False).any(), axis=1)
        else:
            for col in df:
                mask |= df[col].astype(str).str.contains(query, case=False, regex=False, na=False)
        df = df.loc[mask]
    names = list(df.columns)
    if columns.strip():
        terms = [t.strip().lower() for t in columns.split(",") if t.strip()]
        names = [c for c in names if any(t in c.lower() for t in terms)]
        if len(df.columns) and df.columns[0] not in names:
            names.insert(0, df.columns[0])
    return df.loc[:, names[:limit]]


def clear_cache():
    for fn in (_describe, _metadata, _embedding, _expression, _table, _notebook):
        fn.cache_clear()
