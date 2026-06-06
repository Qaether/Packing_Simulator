from __future__ import annotations


def classify_compression_trace(summary_df, pressure_quantile: float = 0.8):
    rows = []
    if summary_df.empty:
        return rows
    p_cut = summary_df["pressure"].quantile(pressure_quantile)
    for _, row in summary_df.iterrows():
        label = "flowing"
        if row["pressure"] >= p_cut and row["msd_step"] < max(summary_df["msd_step"].median(), 1.0e-12):
            label = "jammed_candidate"
        elif row["energy"] > summary_df["energy"].median() and row["msd_step"] < summary_df["msd_step"].median():
            label = "frustrated_candidate"
        rows.append({"phi": row["phi"], "label": label})
    return rows
