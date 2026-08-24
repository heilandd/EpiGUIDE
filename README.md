# EpiGuide NeuroScore

Predict a binary **NeuroScore** (neural-high vs neural-low) from **intraoperative nanopore
sequencing** of a glioma, in seconds, from a single BAM. The NeuroScore is a methylation-derived
marker of neural / non-tumour tissue content; EpiGuide reconstructs it from the sparse per-CpG
methylation calls available during surgery, using a sparse graph neural network
(`SparseForcedEdgeGNN`) over a fixed 4921-CpG similarity graph with masked pooling.

This repository contains the **trained model** (`model/`) and a **minimal BAM→NeuroScore**
pipeline. It accompanies the EpiGuide manuscript.

> Research use only. Not a medical device.

## Install

```bash
git clone <this-repo> && cd Repository
python -m pip install -r requirements.txt
# for BAM input, also install modkit (external): https://github.com/nanoporetech/modkit
```

`torch 2.2` requires **numpy < 2** (already pinned in `requirements.txt`).

## Usage

**From an aligned BAM** (ONT 5mC MM/ML tags, GRCh38/hg38) — end to end:

```bash
python -m neuralscore --bam sample.bam --sample-id PAT001 --out result.csv
```

Internally this runs `modkit pileup` restricted to the model CpG sites
(`model/model_cpgs_hg38.bed`), builds per-CpG calls, and predicts.

**From precomputed per-CpG calls** (no modkit needed):

```bash
python -m neuralscore --calls example/example_calls.tsv --sample-id PAT001 --out result.csv
```

**From Python:**

```python
from neuralscore import load_bundle, predict_from_calls
bundle = load_bundle("model")
print(predict_from_calls(bundle, "example/example_calls.tsv"))
```

### Input: per-CpG calls table
Tab/comma-separated, one row per CpG; duplicate CpGs are summed:

| column | meaning |
|---|---|
| `ID_REF` | Illumina CpG probe id (e.g. `cg23651812`) |
| `methylation_calls` | methylated read count at that CpG |
| `unmethylation_calls` | unmethylated read count |
| `total_calls` | total reads (optional; else meth+unmeth) |

### Output (single-row CSV / dict)
```json
{
  "prob_neural_high": 0.80, "pred_label": "neural_high",
  "confidence": 0.80, "confidence_category": "intermediate_confidence",
  "threshold": 0.5, "n_overlap": 1264, "n_observed_overlap": 1264,
  "observed_fraction_model": 0.26, "mean_coverage_observed": 1.0
}
```
`prob_neural_high` = P(neural-high); the class is neural-high if ≥ `threshold` (0·5).
`confidence` = max(p, 1−p); categories high ≥0·85, intermediate ≥0·65, else low.

## Model
- `SparseForcedEdgeGNN` (2 message-passing layers, hidden 16); 4921 CpG nodes, 98 420 directed
  edges with frozen, task-supervised weights. Node features: methylation call, observed mask,
  coverage. Masked mean pooling over observed nodes + mean over all nodes → MLP → logit.
- Developed on external, dense EPIC methylation-array data with simulated sparse low-coverage
  masking; applied unchanged to intraoperative data. See the manuscript appendix for the full spec.
- Bundle (`model/`): `model_state.pt`, `model_cpgs.csv`, `edge_index.npy`, `edge_attr.npy`,
  `bundle_metadata.json`, `model_cpgs_hg38.bed` (hg38 coordinates of the panel-covered CpGs;
  the rCNS2 intraoperative panel covers ~2600 of the 4921 graph CpGs — the rest are masked, the
  designed low-coverage operating regime).

## Repository layout
```
neuralscore/
  model.py         SparseForcedEdgeGNN + FixedEdgeWeightConv
  bundle.py        load_bundle()
  predict.py       calls_to_data(), predict_from_calls()
  bam_to_calls.py  BAM -> per-CpG calls via modkit
  __main__.py      CLI
model/             trained weights + CpG graph + hg38 bed
example/           example_calls.tsv
```

## Reproducibility
The bundled example reproduces the reference prediction exactly
(`prob_neural_high = 0.7973510`).

## Citation
[Author list]. Intraoperative detection of a methylation-derived NeuroScore from nanopore
sequencing in glioma: the EpiGuide study. *The Lancet* (2026). [DOI].

## License
MIT (see `LICENSE`).
