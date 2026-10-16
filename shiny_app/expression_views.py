"""Build aligned multi-gene tables for plots, summaries, and exports."""
import numpy as np
import pandas as pd


def expression_frame(values, obs, group):
    frames = []
    for gene, vector in values.items():
        frame = pd.DataFrame({"cell_id": obs.cell_id, "gene": gene,
                              "expression": np.asarray(vector)[obs.index]}, index=obs.index)
        frame["group"] = (obs[group].astype("string").fillna("(missing)")
                          if group in obs else "All cells")
        frames.append(frame)
    return pd.concat(frames) if frames else pd.DataFrame()


def expression_summary(frame):
    if frame.empty:
        return pd.DataFrame()
    groups = ["gene", "group"]
    summary = frame.groupby(groups, observed=True).expression.agg(
        cells="size", mean="mean", median="median", maximum="max")
    summary["percent_positive"] = (
        frame.assign(positive=frame.expression > 0)
        .groupby(groups, observed=True).positive.mean() * 100
    )
    return summary.reset_index().round(4)
