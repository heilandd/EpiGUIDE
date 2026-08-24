"""Turn per-CpG intraoperative methylation calls into a NeuroScore prediction.

Input: a table with columns ID_REF (Illumina cg probe id), methylation_calls,
unmethylation_calls, and (optional) total_calls. Duplicate CpGs are summed.
The 3 node features are: methylation call (meth/total), observed mask, and
coverage (min(total, max_coverage)/max_coverage). Unobserved CpGs are masked.
"""
from pathlib import Path
import numpy as np, pandas as pd, torch
from torch_geometric.data import Data


def calls_to_data(path_or_df, bundle, max_coverage=3,
                  cpg_col="ID_REF", meth_col="methylation_calls",
                  unmeth_col="unmethylation_calls", total_col="total_calls"):
    calls = pd.read_csv(path_or_df, sep=None, engine="python") if isinstance(path_or_df, (str, Path)) else path_or_df.copy()
    calls[cpg_col] = calls[cpg_col].astype(str)
    calls[meth_col] = pd.to_numeric(calls[meth_col], errors="coerce").fillna(0)
    calls[unmeth_col] = pd.to_numeric(calls[unmeth_col], errors="coerce").fillna(0)
    if total_col not in calls.columns:
        calls[total_col] = calls[meth_col] + calls[unmeth_col]
    else:
        calls[total_col] = pd.to_numeric(calls[total_col], errors="coerce").fillna(calls[meth_col] + calls[unmeth_col])
    calls = calls.groupby(cpg_col).agg({meth_col: "sum", unmeth_col: "sum", total_col: "sum"}).reset_index()

    model_cpgs = pd.Index(bundle["model_cpgs"])
    cmap = calls.set_index(cpg_col)
    methyl_call = np.zeros(len(model_cpgs), dtype=np.float32)
    observed_mask = np.zeros(len(model_cpgs), dtype=np.float32)
    coverage_norm = np.zeros(len(model_cpgs), dtype=np.float32)

    overlap = model_cpgs.intersection(cmap.index)
    idx = model_cpgs.get_indexer(overlap)
    total = cmap.loc[overlap, total_col].values.astype(np.float32)
    meth = cmap.loc[overlap, meth_col].values.astype(np.float32)
    obs = total > 0
    methyl_call[idx[obs]] = meth[obs] / total[obs]
    observed_mask[idx[obs]] = 1.0
    coverage_norm[idx[obs]] = np.minimum(total[obs], max_coverage) / max_coverage

    x = np.stack([methyl_call, observed_mask, coverage_norm], axis=1).astype(np.float32)
    data = Data(x=torch.from_numpy(x),
                edge_index=torch.as_tensor(np.asarray(bundle["edge_index"]), dtype=torch.long),
                edge_attr=torch.as_tensor(np.asarray(bundle["edge_attr"], dtype=np.float32)))
    data.batch = torch.zeros(data.x.shape[0], dtype=torch.long)
    summary = {"n_cpg_model": len(model_cpgs), "n_overlap": int(len(overlap)),
               "n_observed_overlap": int(obs.sum()), "observed_fraction_model": float(observed_mask.mean()),
               "mean_coverage_observed": float(total[obs].mean()) if obs.sum() else 0.0}
    return data, summary


def _confidence_category(confidence, bundle):
    th = bundle["confidence_thresholds"]
    if confidence >= th.get("high", 0.85): return "high_confidence"
    if confidence >= th.get("intermediate", 0.65): return "intermediate_confidence"
    return "low_confidence"


def predict_from_calls(bundle, path_or_df, max_coverage=3):
    data, summary = calls_to_data(path_or_df, bundle, max_coverage=max_coverage)
    data = data.to(bundle["device"])
    with torch.no_grad():
        prob_high = torch.sigmoid(bundle["model"](data)).item()
    threshold = float(bundle["threshold"])
    pred = int(prob_high >= threshold)
    confidence = max(prob_high, 1 - prob_high)
    return {"prob_neural_high": prob_high, "pred_binary": pred,
            "pred_label": "neural_high" if pred else "neural_low",
            "confidence": confidence, "confidence_category": _confidence_category(confidence, bundle),
            "threshold": threshold, **summary}
