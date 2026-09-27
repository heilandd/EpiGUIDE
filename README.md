<div align="center">

<img src="assets/logo.png" alt="EpiGuide" width="720"/>

<h1>EpiGuide&nbsp;·&nbsp;NeuralScore</h1>

<b>Intraoperative methylation NeuroScore from nanopore sequencing — in seconds, from one BAM.</b>

<br/><br/>

<a href="LICENSE"><img alt="license" src="https://img.shields.io/badge/license-MIT-00468B"></a>
<img alt="python" src="https://img.shields.io/badge/python-3.10%2B-00468B">
<img alt="pytorch" src="https://img.shields.io/badge/PyTorch-2.2-EE4C2C">
<img alt="model" src="https://img.shields.io/badge/model-SparseForcedEdgeGNN-925E9F">
<img alt="graph" src="https://img.shields.io/badge/CpG%20graph-4921%20nodes%20·%2098k%20edges-0099B4">
<img alt="use" src="https://img.shields.io/badge/research%20use%20only-AD002A">

</div>

<br/>

The **NeuralScore** is a methylation-derived marker of neural / non-tumour tissue content.
EpiGuide reconstructs its binary form — **neural-high** vs **neural-low** — from the sparse,
low-coverage per-CpG methylation calls available *during surgery*, using a sparse graph neural
network (`SparseForcedEdgeGNN`) over a fixed 4921-CpG similarity graph with masked pooling.
This repository contains the **trained model** (`model/`) and a **minimal BAM → NeuroScore**
pipeline accompanying the EpiGuide manuscript.

```text
   BAM (ONT 5mC, hg38)  ──modkit──▶  per-CpG calls  ──GNN──▶  NeuralScore + confidence
```

<div align="center"><sub>⚕️ Research use only — not a medical device.</sub></div>

---

## Quick start

```bash
git clone <this-repo> && cd Repository
python -m pip install -r requirements.txt          # torch, torch-geometric, numpy<2, pandas
# for BAM input also install the TESTED modkit (0.3.0 — see VERSIONS.md):
conda install -c bioconda ont-modkit=0.3.0         # https://github.com/nanoporetech/modkit
```

```bash
# end to end, from an aligned nanopore BAM (ONT 5mC MM/ML tags, GRCh38/hg38):
python -m neuralscore --bam sample.bam --sample-id PAT001 --out result.csv

# or from precomputed per-CpG calls (no modkit needed):
python -m neuralscore --calls example/example_calls.tsv --out result.csv
```

```python
from neuralscore import load_bundle, predict_from_calls
bundle = load_bundle("model")
print(predict_from_calls(bundle, "example/example_calls.tsv"))
```

> `torch 2.2` requires **numpy < 2** (already pinned in `requirements.txt`). The BAM path runs
> `modkit pileup` restricted to the model CpG sites (`model/model_cpgs_hg38.bed`), then predicts.
> **modkit versions matter** — tested with **modkit 0.3.0**; newer releases changed the
> `pileup --include-bed` CLI. Exact versions, the tested command, and an automatic
> compatibility fallback are documented in [`VERSIONS.md`](VERSIONS.md).

## Input — per-CpG calls table

Tab- or comma-separated, one row per CpG (duplicate CpGs are summed):

| column | meaning |
|---|---|
| `ID_REF` | Illumina CpG probe id (e.g. `cg23651812`) |
| `methylation_calls` | methylated read count at that CpG |
| `unmethylation_calls` | unmethylated read count |
| `total_calls` | total reads *(optional; else meth + unmeth)* |

## Output

```json
{
  "prob_neural_high": 0.80, "pred_label": "neural_high",
  "confidence": 0.80, "confidence_category": "intermediate_confidence",
  "threshold": 0.5, "n_overlap": 1264, "n_observed_overlap": 1264,
  "observed_fraction_model": 0.26, "mean_coverage_observed": 1.0
}
```

`prob_neural_high` = P(neural-high); the class is neural-high if ≥ `threshold` (0·5).
`confidence` = max(p, 1−p); categories: **high** ≥ 0·85, **intermediate** ≥ 0·65, else **low**.

## Model

- **`SparseForcedEdgeGNN`** — 2 message-passing layers (hidden 16) over **4921 CpG nodes** and
  **98 420** directed edges with frozen, task-supervised weights. Node features: methylation
  call, observed mask, coverage. Masked mean pooling over observed nodes + mean over all nodes
  → MLP → logit.
- Developed on external, dense EPIC methylation-array data with simulated sparse low-coverage
  masking, then applied **unchanged** to intraoperative data (full spec in the manuscript appendix).
- **Bundle** (`model/`): `model_state.pt`, `model_cpgs.csv`, `edge_index.npy`, `edge_attr.npy`,
  `bundle_metadata.json`, and `model_cpgs_hg38.bed` — hg38 coordinates of the panel-covered CpGs.
  The rCNS2 intraoperative panel covers ≈ 2600 of the 4921 graph CpGs; the rest are masked — the
  designed low-coverage operating regime.

## Layout

```text
neuralscore/
  model.py         SparseForcedEdgeGNN + FixedEdgeWeightConv
  bundle.py        load_bundle()
  predict.py       calls_to_data(), predict_from_calls()
  bam_to_calls.py  BAM → per-CpG calls via modkit
  __main__.py      CLI (--bam | --calls)
model/             trained weights + CpG graph + hg38 bed
example/           example_calls.tsv
```

## Reproducibility

The bundled example reproduces the reference prediction exactly — `prob_neural_high = 0.7973510`.

## Citation

Hope comming soon!

## ⚖️ License

Code released under the **MIT License** (see [`LICENSE`](LICENSE)). Research use only.

<div align="center"><br/><sub>EpiGuide&nbsp;·&nbsp;Universitätsklinikum Erlangen (FAU)&nbsp;·&nbsp;UKE Hamburg</sub></div>
