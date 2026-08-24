"""Load the released model bundle (weights, CpG graph, thresholds)."""
import json
from pathlib import Path
import numpy as np, pandas as pd, torch
from .model import SparseForcedEdgeGNN


def load_bundle(bundle_dir, device=None):
    bundle_dir = Path(bundle_dir)
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    meta = json.load(open(bundle_dir / "bundle_metadata.json"))

    model = SparseForcedEdgeGNN(**meta["model_config"])
    model.load_state_dict(torch.load(bundle_dir / "model_state.pt", map_location=device))
    model.to(device).eval()

    return {
        "model": model,
        "edge_index": np.load(bundle_dir / "edge_index.npy"),
        "edge_attr": np.load(bundle_dir / "edge_attr.npy"),
        "model_cpgs": pd.read_csv(bundle_dir / "model_cpgs.csv")["CpGs"].astype(str).tolist(),
        "threshold": float(meta["threshold"]),
        "confidence_thresholds": meta["confidence_thresholds"],
        "metadata": meta,
        "device": device,
    }
