# EDNA single-cell explorer

A local, read-only **Python Shiny** app connected to the project folder. The app discovers saved files recursively and refreshes when files are added or changed. The existing notebooks and scientific data are left unchanged.

## Run

Double-click `shiny_app/launch.command`, or run from the project folder:

```bash
shiny_app/.venv/bin/python -m shiny run --host 127.0.0.1 --port 8765 shiny_app/app.py
```

Open http://127.0.0.1:8765. Keep the terminal running; press Control-C to stop. If that port is already in use, use `--port 8766`, or `EDNA_PORT=8766 ./shiny_app/launch.command`.

The environment installed on this computer uses its existing scientific Python packages through `--system-site-packages`; Shiny and plotting packages are installed locally in `.venv`. To create a self-contained environment on another computer:

```bash
python3 -m venv shiny_app/.venv
shiny_app/.venv/bin/python -m pip install -r shiny_app/requirements.txt
```

Python 3.11–3.13 is recommended. Set `EDNA_PROJECT_ROOT` before launch to browse another project folder. The app binds to loopback and has no authentication; use an authenticated Shiny hosting setup if remote access is needed.

## Views

- **Overview:** active cell and gene totals, UMAP (or the first available embedding), condition counts, and recent saved outputs.
- **Tables:** CSV/TSV, Excel sheets, Parquet, AnnData cell/gene metadata, saved `uns` tables/arrays, embeddings, 10x feature/barcode metadata, matrix inventories, and small expression previews. Search rows or column names, sort/filter the grid, export the search result, or download the source.
- **Cells & QC:** UMAP/PCA colored by metadata, condition and singlet filters, numeric QC distributions, composition counts, and filtered metadata export.
- **Gene expression:** select up to six genes to compare separate UMAP panels with a shared expression color scale and the same cells in every panel. Choose `X`, `layers/*`, or `raw/X`; distributions, summaries, and CSV exports include each selected gene separately. Values are used as stored; the app does not normalize or combine genes. Percentage positive means values greater than zero, which may differ from detection rates for transformed or centered layers.
- **Figures:** saved raster figures and downloadable PDFs/vectors. The notebook-output viewer, project-file browser, and recent-file list are not shown.
- **Colors:** 24 distinct dark categorical colors; numeric expression uses a high-contrast gray → blue → purple → red → gold scale. Multi-gene maps share the same scale.

## Current clustering data

The default dataset is `results_intermediate/clustered_filtered.h5ad`: **20,403 cells × 14,331 genes**, with **15 clusters (0–14)**. Clusters 15 and 16 (20 cells total) have been removed. The original `decontaminated/pp_and_umap_data.h5ad` is preserved.

This initial app dataset was recreated from that preprocessed file using the notebook's igraph Leiden settings (resolution 0.5, default random seed 0), with labels ranked by cluster size. Every cluster-by-condition count matched the saved output of notebook execution [42] before removal. It is a recomputed dataset, not a capture of the live notebook session. The existing embeddings are retained, not recomputed after filtering. Provenance is stored in `uns/shiny_export`.

The Overview colors UMAP by cluster when available and shows a cluster-by-condition count table. Sidebar cluster, condition, and singlet filters apply to cell views and their exports. Cells & QC and Gene expression default to comparisons by cluster. `results_intermediate/annotation_guess.csv` remains an annotation hypothesis table; its guesses are not automatically assigned as confirmed cell types.

## Update from the notebook

Run the final **Save the current filtered analysis for the Shiny explorer** cell in `clustering_and_annotation.ipynb` after your clustering, filtering, and any annotation work. It saves the current `adata` to `results_intermediate/clustered_filtered.h5ad` using a temporary file followed by rename, and exports `cell_counts_by_cluster_condition.csv`. The export checks that clusters 15 and 16 have been removed. Any marker tables stored in `adata.uns` and cell-type columns stored in `adata.obs` are included in the saved object.

The app watches saved files every 20 seconds. Click **Refresh project** to load changes immediately. In an already open session, select `clustered_filtered.h5ad` under **Active dataset**; fresh sessions select it by default. The app does not read unsaved notebook memory or execute notebook code.

## Limits and data interpretation

- Cell filters affect Overview, Cells & QC, and Gene expression. Tables use their selected source independently.
- Table previews load up to **100,000 rows** for CSV/TSV/Excel/Parquet and 10x barcode lists. At most **60 columns** are displayed; search column names to choose others. AnnData metadata is read in full for cell filtering. Original-file downloads contain the full source.
- **Export this view** includes the row-search and column-search result; additional grid header filters and grid sorting only affect the on-screen grid. Cell and expression downloads include all cells selected by the sidebar filters.
- Expression-matrix previews show the first **50 cells × 20 genes**. They do not represent the full matrix. Large matrices are never densified in full; expression plots read one gene at a time in sparse chunks.
- Overview and gene plots show up to **30,000 cells**, with a deterministic random sample above that. Cells & QC lets you choose its plotting limit. Counts, summary statistics, and exports always use the full current cell selection.
- A file being written may temporarily fail to read. Refresh after its save completes. Closed HDF5 handles and caches keyed by file modification time prevent stale handle reuse.
- The app does not infer biology, modify annotations, rerun preprocessing, or execute notebook code. `scDblFinder.cluster` is a preprocessing field, not a substitute for final Leiden clusters.

## Verify

```bash
shiny_app/.venv/bin/python -m unittest discover -s shiny_app/tests -v
```

The tests compare matrix readers against known dense values, exercise CSR/CSC/dense storage, raw gene axes, file boundaries, cache invalidation, notebook snapshots, and literal table searches. See [Shiny's table API](https://shiny.posit.co/py/api/core/render.data_frame.html) for the underlying grid behavior.

## Share with collaborators

The current `http://127.0.0.1:8765` address is local to your computer; sending that link does not give others access.

For a hosted URL, deploy the Python app to **Posit Connect Cloud**, or your institution's **Posit Connect** server if available. Connect Cloud supports Python Shiny deployment using Posit Publisher or the rsconnect-python CLI. Choose the access setting that fits your audience; private sharing and named-user access depend on the hosting plan.

Prepare a separate deployment folder containing only the app and intended datasets:

```text
edna-share/
  requirements.txt                  # copy shiny_app/requirements.txt
  shiny_app/
    app.py
    data_access.py
    expression_views.py
    www/style.css
  results_intermediate/
    clustered_filtered.h5ad
    annotation_guess.csv            # optional
    cell_counts_by_cluster_condition.csv  # optional
```

Use `shiny_app/app.py` as the entry point, Python 3.11–3.13, and the top-level `requirements.txt`. This structure keeps the default project root correct. Include any figures you want in a separate figures folder. Do not include `.venv`, notebooks, raw downloads, or the rest of the working project. Tables and their source downloads remain available, so include only files intended for viewers. Removing navigation tabs is not an access-control mechanism.

The current AnnData file is about **827 MiB**. Include it directly in the deployment bundle with Posit Publisher; it exceeds GitHub's ordinary per-file limit. Check your host's bundle-size, storage, and memory allowances before uploading. Hosted data is a snapshot: local notebook saves do not update the hosted copy; redeploy the updated data.

Alternatively, share this curated folder privately with collaborators who can install Python. From `edna-share`, they can create a virtual environment, install `requirements.txt`, and run `python -m shiny run shiny_app/app.py` locally.

References:
- https://docs.posit.co/connect-cloud/how-to/python/shiny-python.html
- https://docs.posit.co/connect-cloud/user/share/index.html
- https://shiny.posit.co/py/docs/deploy-on-prem.html
