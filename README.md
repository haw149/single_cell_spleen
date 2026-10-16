# Single-cell EDNA analysis and Shiny explorer

## Deployment settings

- Repository: `haw149/single_cell_spleen`
- Branch: `main`
- Shiny primary file: `shiny_app/app.py`
- Python dependencies: root `requirements.txt`
- Main analysis notebook: `clustering_and_annotation.ipynb`

## Notebooks

Run preprocessing in order: `preprocessing-1.ipynb`, then `preprocessing-2.ipynb`, followed by `clustering_and_annotation.ipynb`. Notebooks include their saved outputs and currently use local data paths; adapt those paths for another machine. The Shiny requirements cover the explorer, not all notebook analysis dependencies.

## Run the explorer

From this repository root, create a Python 3.11–3.13 environment, install `requirements.txt`, and run:

```bash
python -m pip install -r requirements.txt
python -m shiny run --host 127.0.0.1 --port 8765 shiny_app/app.py
```

Open http://127.0.0.1:8765.

## Required data

The large scientific datasets are not stored in GitHub. Supply `results_intermediate/clustered_filtered.h5ad` in the deployment bundle or on the hosting server. The current local export contains 20,403 cells and 14,331 genes, with clusters 15 and 16 removed. It is about 827 MiB. Optionally supply annotation/count CSVs in `results_intermediate/` and saved figures.

The source code alone will launch an empty explorer until data is provided. For Posit hosting, select `shiny_app/app.py` as the entry point and include the dataset separately. Do not upload the entire working directory. See [app documentation](shiny_app/README.md) for details.

## Verify

```bash
python -m unittest discover -s shiny_app/tests -v
```
