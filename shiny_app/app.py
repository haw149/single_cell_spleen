"""EDNA — a read-only Shiny browser for the single-cell project."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from shiny import App, reactive, render, ui
from shinywidgets import output_widget, render_plotly

sys.path.insert(0, str(Path(__file__).parent))
import data_access as data
import expression_views

HERE = Path(__file__).parent
PALETTE = px.colors.qualitative.Dark24
EXPRESSION_SCALE = [
    [0.0, "#d1d5db"], [0.15, "#2563eb"], [0.45, "#6d28d9"],
    [0.75, "#be123c"], [1.0, "#ffb000"],
]


def note(text, error=False):
    return ui.div(text, class_="notice error-box" if error else "notice")


def styled(fig, height=390):
    fig.update_layout(template="plotly_white", font=dict(family="Inter, -apple-system, sans-serif", size=11, color="#243244"),
                      margin=dict(l=30, r=15, t=25, b=35), height=height, paper_bgcolor="white", plot_bgcolor="white",
                      legend=dict(title_text="", orientation="h", y=-0.15), colorway=PALETTE)
    fig.update_xaxes(gridcolor="#e2e8f0", zeroline=False)
    fig.update_yaxes(gridcolor="#e2e8f0", zeroline=False)
    return fig


def empty_plot(message):
    fig = styled(go.Figure())
    fig.add_annotation(text=message, x=.5, y=.5, xref="paper", yref="paper", showarrow=False)
    fig.update_xaxes(visible=False); fig.update_yaxes(visible=False)
    return fig


def grid(df, height="470px"):
    df = df.copy()
    for col in df:
        if isinstance(df[col].dtype, pd.CategoricalDtype):
            df[col] = df[col].astype("string").fillna("(missing)")
    return render.DataGrid(df, width="100%", height=height, filters=True, selection_mode="rows")


app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.div(ui.div(ui.span("◌", class_="brand-mark"), "EDNA", class_="brand"),
               ui.div("SINGLE-CELL WORKSPACE", class_="brand-sub")),
        ui.div(ui.div("Your project", class_="eyebrow"), ui.p("An open view of your analysis.", class_="muted")),
        ui.input_select("dataset", "ACTIVE DATASET", choices={}),
        ui.output_ui("dataset_path"),
        ui.input_selectize("conditions", "Conditions · cell views", choices=[], multiple=True,
                           options={"placeholder": "All conditions"}),
        ui.input_selectize("clusters", "Clusters · cell views", choices=[], multiple=True,
                           options={"placeholder": "All clusters"}),
        ui.input_switch("singlets", "Singlets only", value=False),
        ui.input_action_button("refresh", "↻  Refresh project", class_="btn-primary w-100"),
        ui.input_switch("auto_refresh", "Watch for saved outputs", value=True),
        ui.div(ui.tags.strong("Connected to your files"), ui.p("New tables, datasets and figures appear as you save them. The watcher checks every 20 seconds."),
               ui.span("Cell filters apply to Overview, Cells & QC and Gene expression. Tables show their own source."), class_="sidebar-note"),
        width=255, open="desktop",
    ),
    ui.tags.head(ui.include_css(HERE / "www" / "style.css")),
    ui.div(ui.div(ui.div("EXPERIMENT / EDNA", class_="eyebrow"), ui.h1("Your data, in view."),
                  ui.p("Explore the cells, tables and results behind your single-cell analysis.")),
           ui.span("Local project", class_="live-badge"), class_="hero"),
    ui.output_ui("stats"),
    ui.navset_pill(
        ui.nav_panel("Overview",
            ui.output_ui("dataset_notice"),
            ui.layout_columns(
                ui.card(ui.card_header("Cell landscape", ui.span("  / UMAP", class_="muted")), output_widget("overview_umap")),
                ui.card(ui.card_header("Cells by condition"), output_widget("condition_plot")), col_widths=[7, 5]),
            ui.card(ui.card_header("Inside the active dataset"), ui.output_ui("dataset_summary")),
            ui.card(ui.card_header("Cells by cluster and condition"), ui.output_data_frame("cluster_counts")),
        ),
        ui.nav_panel("Tables",
            ui.card(ui.card_header("Table browser"),
                ui.div(ui.input_select("table_file", "Source file", {}), ui.input_select("table_key", "Table or matrix", {}), class_="table-toolbar"),
                ui.div(ui.input_text("row_search", "Search all loaded columns", placeholder="A gene, barcode, cluster, or value…"),
                       ui.input_text("column_search", "Show columns containing", placeholder="All columns · or e.g. gene, counts, mt"),
                       ui.download_button("download_table", "↓ Export this view"), class_="search-toolbar"),
                ui.output_ui("table_info"), ui.output_data_frame("table_grid"),
                ui.download_button("download_table_source", "Download original source file")),
        ),
        ui.nav_panel("Cells & QC",
            ui.card(ui.card_header("Explore individual cells"),
                ui.div(ui.input_select("embedding", "Embedding", {}), ui.input_select("cell_color", "Color by", {}),
                       ui.input_slider("plot_limit", "Maximum points shown", min=1000, max=50000, value=30000, step=1000), class_="plot-controls"),
                output_widget("cell_plot"), ui.output_ui("cell_plot_note")),
            ui.layout_columns(
                ui.card(ui.card_header("Quality-control distribution"), ui.input_select("qc_metric", "Metric", {}), output_widget("qc_plot")),
                ui.card(ui.card_header("Cell composition"), ui.input_select("composition_group", "Group cells by", {}), output_widget("composition_plot")), col_widths=[6, 6]),
            ui.card(ui.card_header("Cell metadata · current cell filters"), ui.download_button("download_cells", "↓ Export filtered cells"),
                    ui.output_data_frame("cell_grid")),
        ),
        ui.nav_panel("Gene expression",
            ui.card(ui.card_header("Gene expression across cells"),
                ui.div(ui.input_selectize("gene", "Genes · up to 6", {}, multiple=True,
                                         options={"maxItems": 6, "placeholder": "Search and select genes"}), ui.input_select("layer", "Expression values", {}),
                       ui.input_select("expression_group", "Compare by", {}), class_="plot-controls"),
                note("Values are shown as stored in the selected matrix. X is the saved analysis matrix; layers/counts contains counts when available. No normalization is applied here."),
                ui.output_ui("expression_status"),
                output_widget("gene_plot"), output_widget("gene_distribution")),
            ui.card(ui.card_header("Expression summary · current cell filters"), ui.output_data_frame("expression_summary"),
                    ui.download_button("download_expression", "↓ Export expression by cell")),
        ),
        ui.nav_panel("Figures",
            ui.card(ui.input_select("figure_file", "Figure", {}),
                    ui.output_ui("figure_status"), ui.output_image("figure_preview", height="auto"),
                    ui.download_button("download_figure", "↓ Download figure")),
        ), id="main_tab", selected="Overview",
    ),
    ui.div(ui.span("EDNA / SINGLE-CELL EXPLORER"), ui.span("Read-only · backed by your project files"), class_="footer"),
    title="EDNA · Single-cell explorer", window_title="EDNA · Single-cell explorer", fillable=False,
)


def server(input, output, session):
    catalog = reactive.value(data.discover_files())
    revision = reactive.value(0)

    @reactive.effect
    def watch_files():
        input.refresh()
        if input.auto_refresh(): reactive.invalidate_later(20)
        new = data.discover_files()
        with reactive.isolate():
            if not new.equals(catalog.get()):
                catalog.set(new)
                revision.set(revision.get() + 1)

    @reactive.effect
    @reactive.event(input.refresh, ignore_init=True)
    def force_refresh():
        data.clear_cache()
        revision.set(revision.get() + 1)

    def update_choice(id, choices, preferred=None):
        with reactive.isolate(): current = getattr(input, id)()
        selected = current if current in choices else preferred if preferred in choices else next(iter(choices), "")
        ui.update_select(id, choices=choices, selected=selected)

    @reactive.effect
    def update_catalog():
        files = catalog.get()
        choices = {p: Path(p).name for p in files.loc[files.kind == "AnnData", "path"]}
        preferred = ("results_intermediate/clustered_filtered.h5ad" if "results_intermediate/clustered_filtered.h5ad" in choices
                     else "decontaminated/pp_and_umap_data.h5ad")
        update_choice("dataset", choices, preferred)
        for id, kinds in [("table_file", ["AnnData", "HDF5", "Table"]), ("figure_file", ["Figure"])]:
            choices = {p: p for p in files.loc[files.kind.isin(kinds), "path"]}
            preferred = "results_intermediate/annotation_guess.csv" if id == "table_file" else None
            update_choice(id, choices, preferred)

    @reactive.calc
    def metadata():
        revision.get()
        if not input.dataset(): return {"error": "No AnnData file found. Save an .h5ad file in the project folder, then refresh."}
        try: return data.describe_h5ad(input.dataset())
        except Exception as exc: return {"error": f"This dataset could not be read: {exc}. If it is being written, try Refresh when the save finishes."}

    @reactive.calc
    def all_cells():
        revision.get()
        if "error" in metadata(): return pd.DataFrame()
        try: return data.load_obs(input.dataset())
        except Exception: return pd.DataFrame()

    @reactive.effect
    def update_dataset_controls():
        meta, obs = metadata(), all_cells()
        categorical = {c: c for c in obs if c != "cell_id" and (not pd.api.types.is_numeric_dtype(obs[c]) or obs[c].nunique() <= 40)}
        numeric = {c: c for c in obs if pd.api.types.is_numeric_dtype(obs[c]) and not pd.api.types.is_bool_dtype(obs[c])}
        update_choice("embedding", {k: k.removeprefix("X_").upper() for k in meta.get("embeddings", [])}, "X_umap")
        group = "cluster" if "cluster" in obs else "condition"
        update_choice("cell_color", {**categorical, **numeric}, group)
        update_choice("composition_group", categorical, group)
        update_choice("expression_group", categorical, group)
        update_choice("qc_metric", numeric, "total_counts")
        update_choice("layer", {k: k for k in meta.get("layers", [])}, "X")
        choices = sorted(obs["condition"].dropna().astype(str).unique()) if "condition" in obs else []
        with reactive.isolate(): selected = [c for c in (input.conditions() or []) if c in choices]
        ui.update_selectize("conditions", choices=choices, selected=selected)
        choices = sorted(obs["cluster"].dropna().astype(str).unique(),
                         key=lambda x: (0, int(x)) if x.isdigit() else (1, x)) if "cluster" in obs else []
        with reactive.isolate(): selected = [c for c in (input.clusters() or []) if c in choices]
        ui.update_selectize("clusters", choices=choices, selected=selected)

    @reactive.effect
    def update_genes():
        revision.get()
        if not input.dataset() or not input.layer(): return
        try:
            genes = data.gene_names(input.dataset(), input.layer())
            with reactive.isolate(): current = input.gene()
            chosen = [g for g in (current or []) if g in genes]
            if not chosen: chosen = ["Cd3d"] if "Cd3d" in genes else genes[:1]
            ui.update_selectize("gene", choices=genes, selected=chosen, server=True)
        except Exception: ui.update_selectize("gene", choices=[], selected=[])

    @reactive.calc
    def cells():
        obs = all_cells()
        if "condition" in obs and input.conditions(): obs = obs.loc[obs.condition.astype(str).isin(input.conditions())]
        if "cluster" in obs and input.clusters(): obs = obs.loc[obs["cluster"].astype(str).isin(input.clusters())]
        if input.singlets() and "scDblFinder.class" in obs:
            obs = obs.loc[obs["scDblFinder.class"].astype(str).str.lower().eq("singlet")]
        return obs

    @render.ui
    def dataset_path():
        return ui.div(input.dataset() or "No dataset saved yet", class_="source-path")

    @render.ui
    def stats():
        meta, obs = metadata(), cells()
        items = [("ACTIVE CELLS", f"{len(obs):,}", f"of {meta.get('n_obs', 0):,} saved cells"),
                 ("GENES", f"{meta.get('n_vars', 0):,}", "in the active dataset"),
                 ("CONDITIONS", str(obs.condition.nunique()) if "condition" in obs else "—", "in the current cell selection"),
                 ("CLUSTERS", str(obs["cluster"].nunique()) if "cluster" in obs else "—", "in the current cell selection")]
        return ui.div(*[ui.div(ui.div(label, class_="eyebrow"), ui.div(value, class_="stat-value"), ui.div(caption, class_="stat-label"),
                              class_="stat featured" if i == 0 else "stat") for i,(label,value,caption) in enumerate(items)], class_="stat-grid")

    @render.ui
    def dataset_notice():
        meta = metadata()
        if "error" in meta: return note(meta["error"], True)
        if not any("cluster" == c or "leiden" in c.lower() or "cell_type" in c.lower() for c in meta["obs_columns"]):
            return note("Final clustering and annotation labels are not saved in this dataset yet. Annotation CSVs are available in Tables. Save the annotated .h5ad to make its labels available in cell views.")
        return None

    def scatter(embedding, color, frame, values=None, height=440, limit=30000):
        if not embedding or frame.empty: return empty_plot("No embedding or cells available for this selection.")
        try:
            coords = data.read_embedding(input.dataset(), embedding)
            plot = pd.DataFrame({"x": coords[frame.index, 0], "y": coords[frame.index, 1], "cell_id": frame.cell_id.to_numpy()}, index=frame.index)
            if values is not None:
                plot[color] = np.asarray(values)[frame.index]
            elif color in frame:
                series = frame[color]
                plot[color] = series.to_numpy() if pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series) else series.astype("string").fillna("(missing)").to_numpy()
            else: color = None
            if len(plot) > limit: plot = plot.sample(limit, random_state=42)
            if color and pd.api.types.is_numeric_dtype(plot[color]): plot = plot.sort_values(color)
            fig = px.scatter(plot, x="x", y="y", color=color, hover_name="cell_id", render_mode="webgl",
                             color_discrete_sequence=PALETTE, color_continuous_scale=EXPRESSION_SCALE,
                             labels={"x": embedding.removeprefix("X_").upper()+" 1", "y": embedding.removeprefix("X_").upper()+" 2"})
            fig.update_traces(marker=dict(size=3, opacity=.9))
            return styled(fig, height)
        except Exception as exc: return empty_plot(f"Unable to read embedding: {str(exc)[:180]}")

    @render_plotly
    def overview_umap():
        embeddings = metadata().get("embeddings", [])
        embedding = "X_umap" if "X_umap" in embeddings else embeddings[0] if embeddings else None
        return scatter(embedding, "cluster" if "cluster" in cells() else "condition", cells(), height=370)

    @render.data_frame
    def cluster_counts():
        obs = cells()
        if not {"cluster", "condition"}.issubset(obs.columns):
            return grid(pd.DataFrame({"Status": ["Save cluster and condition labels to see counts here."]}), "180px")
        counts = pd.crosstab(obs["cluster"].astype(str), obs["condition"].astype(str))
        counts = counts.loc[sorted(counts.index, key=lambda x: (0, int(x)) if x.isdigit() else (1, x))]
        counts["Total"] = counts.sum(axis=1)
        return grid(counts.reset_index(), "360px")

    def composition(column, frame, height=370):
        if column not in frame or frame.empty: return empty_plot("No grouping information available.")
        counts = frame[column].astype("string").fillna("(missing)").value_counts().rename_axis("Group").reset_index(name="Cells")
        counts = counts.sort_values("Cells")
        fig = px.bar(counts, x="Cells", y="Group", orientation="h", color="Group", color_discrete_sequence=PALETTE, text="Cells")
        fig.update_traces(textposition="outside", cliponaxis=False)
        styled(fig, height); fig.update_layout(showlegend=False, margin=dict(l=10,r=40,t=25,b=35)); fig.update_yaxes(title=None)
        fig.update_xaxes(range=[0, counts.Cells.max()*1.22])
        return fig

    @render_plotly
    def condition_plot(): return composition("condition", cells())

    @render.ui
    def dataset_summary():
        meta = metadata()
        if "error" in meta: return note(meta["error"], True)
        return ui.TagList(ui.p(f"{len(meta['obs_columns'])} cell metadata fields and {len(meta['var_columns'])} gene metadata fields."),
                         ui.div("EMBEDDINGS", class_="eyebrow"), ui.p(", ".join(meta["embeddings"]) or "None saved"),
                         ui.div("EXPRESSION MATRICES", class_="eyebrow"), ui.p(", ".join(meta["layers"])))

    @reactive.effect
    def update_table_keys():
        revision.get()
        try: choices = data.table_choices(input.table_file()) if input.table_file() else {}
        except Exception: choices = {}
        update_choice("table_key", choices)

    @reactive.calc
    def source_table():
        revision.get()
        try:
            if not input.table_file() or not input.table_key(): return pd.DataFrame(), "Select a source file and a table."
            df = data.load_table(input.table_file(), input.table_key())
            return df, df.attrs.get("note", "")
        except Exception as exc: return pd.DataFrame(), f"Unable to read this table: {exc}. If the file is being saved, refresh when it finishes."

    @reactive.calc
    def displayed_table():
        return data.filter_table(source_table()[0], input.row_search(), input.column_search())

    @render.ui
    def table_info():
        df, message = source_table(); shown = displayed_table()
        return ui.TagList(note(message) if message else None,
            ui.div(f"{len(shown):,} matching rows · {len(shown.columns):,} visible columns / {len(df.columns):,} loaded columns. Up to 60 columns shown; use the column search to find genes or fields. Export includes these matching rows and visible columns; grid header filters only change the on-screen grid.", class_="table-note"))

    @render.data_frame
    def table_grid(): return grid(displayed_table())

    @render.download(filename=lambda: f"{Path(input.table_file() or 'table').stem}_view.csv")
    def download_table(): yield displayed_table().to_csv(index=False)

    @render.download(filename=lambda: Path(input.table_file()).name)
    def download_table_source(): return str(data.safe_path(input.table_file()))

    @render_plotly
    def cell_plot(): return scatter(input.embedding(), input.cell_color(), cells(), limit=input.plot_limit())

    @render.ui
    def cell_plot_note():
        message = f"Showing {min(len(cells()), input.plot_limit()):,} of {len(cells()):,} selected cells. Larger selections use a reproducible random sample for plotting; tables and summaries use every selected cell."
        if input.singlets() and "scDblFinder.class" not in all_cells(): message += " Singlet status is not saved in this dataset; the singlets filter cannot be applied."
        return ui.div(message, class_="table-note")

    @render_plotly
    def qc_plot():
        obs, metric = cells(), input.qc_metric()
        if metric not in obs or obs.empty: return empty_plot("No numeric QC metric available.")
        plot = obs.copy()
        if "condition" in plot: plot["condition"] = plot.condition.astype("string").fillna("(missing)")
        return styled(px.histogram(plot, x=metric, color="condition" if "condition" in plot else None, nbins=60,
                                  barmode="overlay", opacity=.65, color_discrete_sequence=PALETTE))

    @render_plotly
    def composition_plot(): return composition(input.composition_group(), cells())

    @render.data_frame
    def cell_grid(): return grid(cells(), "390px")

    @render.download(filename="filtered_cell_metadata.csv")
    def download_cells(): yield cells().to_csv(index=False)

    @reactive.calc
    def expression():
        revision.get()
        if not input.gene() or not input.layer(): return None, "Select a gene and expression matrix."
        if len(input.gene()) > 6: return None, "Select up to six genes at a time."
        try: return {g: data.read_expression(input.dataset(), g, input.layer()) for g in input.gene()}, ""
        except Exception as exc: return None, f"Unable to read expression: {exc}"

    @reactive.calc
    def expression_frame():
        values, _ = expression()
        if values is None: return pd.DataFrame()
        return expression_views.expression_frame(values, cells(), input.expression_group())

    @render.ui
    def expression_status():
        values, error = expression()
        return note(error, True) if error else ui.div(
            f"{', '.join(input.gene())} · {input.layer()} · one panel per gene, shared expression color scale. "
            f"The same {min(len(cells()), 30000):,} cells appear in every map; all {len(cells()):,} selected cells are included in summaries and exports.", class_="table-note")

    @render_plotly
    def gene_plot():
        values, error = expression()
        if values is None: return empty_plot(error)
        embeddings = metadata().get("embeddings", [])
        embedding = "X_umap" if "X_umap" in embeddings else embeddings[0] if embeddings else None
        frame = expression_frame()
        if not embedding or frame.empty: return empty_plot("No embedding or cells available for this selection.")
        try:
            coords = data.read_embedding(input.dataset(), embedding)
            selected = cells()
            if len(selected) > 30000: selected = selected.sample(30000, random_state=42)
            plot = frame.loc[frame.index.isin(selected.index)].copy()
            plot["x"] = coords[plot.index, 0]
            plot["y"] = coords[plot.index, 1]
            plot = plot.sort_values("expression")
            fig = px.scatter(plot, x="x", y="y", color="expression", facet_col="gene", facet_col_wrap=3,
                             category_orders={"gene": list(values)}, hover_name="cell_id", render_mode="webgl",
                             color_continuous_scale=EXPRESSION_SCALE,
                             labels={"x": embedding.removeprefix("X_").upper()+" 1", "y": embedding.removeprefix("X_").upper()+" 2"})
            fig.update_traces(marker=dict(size=3, opacity=.9))
            return styled(fig, 390 * ((len(values) + 2) // 3))
        except Exception as exc: return empty_plot(f"Unable to read embedding: {str(exc)[:180]}")

    @render_plotly
    def gene_distribution():
        frame = expression_frame()
        if frame.empty: return empty_plot("No expression values for this selection.")
        fig = px.box(frame, x="group", y="expression", color="group", points=False, color_discrete_sequence=PALETTE,
                     facet_col="gene", facet_col_wrap=3, category_orders={"gene": list(input.gene())},
                     labels={"group": input.expression_group() or "Group", "expression": "Expression · " + input.layer()})
        styled(fig, 360 * ((len(input.gene()) + 2) // 3)); fig.update_layout(showlegend=False)
        return fig

    @render.data_frame
    def expression_summary():
        frame = expression_frame()
        if frame.empty: return grid(pd.DataFrame())
        return grid(expression_views.expression_summary(frame), "300px")

    @render.download(filename="selected_genes_expression.csv")
    def download_expression():
        frame = expression_frame().copy(); frame["matrix"] = input.layer()
        yield frame.to_csv(index=False)

    @render.ui
    def figure_status():
        if not input.figure_file(): return note("No saved figures found yet.")
        return ui.div(input.figure_file(), class_="source-path") if Path(input.figure_file()).suffix.lower() in (".png", ".jpg", ".jpeg") else note("Download this vector or PDF figure to open it in your preferred viewer.")

    @render.image
    def figure_preview():
        revision.get()
        if not input.figure_file() or Path(input.figure_file()).suffix.lower() not in (".png", ".jpg", ".jpeg"): return None
        return {"src": str(data.safe_path(input.figure_file())), "alt": Path(input.figure_file()).name, "width": "100%", "class": "figure-preview"}

    @render.download(filename=lambda: Path(input.figure_file()).name)
    def download_figure(): return str(data.safe_path(input.figure_file()))



app = App(app_ui, server)
