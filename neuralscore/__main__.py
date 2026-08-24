"""Command-line: BAM (or precomputed calls) -> NeuroScore.

  python -m neuralscore --bam sample.bam --out result.csv        # end-to-end (needs modkit)
  python -m neuralscore --calls calls.tsv --out result.csv       # from per-CpG calls
"""
import argparse, json, os, sys
import pandas as pd
from . import load_bundle, predict_from_calls
from .bam_to_calls import bam_to_calls

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    ap = argparse.ArgumentParser(description="EpiGuide NeuroScore prediction")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--bam", help="aligned nanopore BAM with 5mC MM/ML tags (hg38)")
    g.add_argument("--calls", help="precomputed per-CpG calls table (ID_REF, methylation_calls, unmethylation_calls, total_calls)")
    ap.add_argument("--model", default=os.path.join(HERE, "model"), help="model bundle directory")
    ap.add_argument("--sample-id", default="sample")
    ap.add_argument("--out", help="output CSV (single-row result)")
    ap.add_argument("--ref", help="reference FASTA for modkit (optional)")
    ap.add_argument("--modkit", default="modkit")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()

    bundle = load_bundle(a.model)
    if a.bam:
        model_bed = os.path.join(a.model, "model_cpgs_hg38.bed")
        calls = bam_to_calls(a.bam, model_bed, modkit=a.modkit, threads=a.threads, ref=a.ref)
        result = predict_from_calls(bundle, calls)
    else:
        result = predict_from_calls(bundle, a.calls)

    result = {"sample_id": a.sample_id, **result}
    print(json.dumps(result, indent=2))
    if a.out:
        pd.DataFrame([result]).to_csv(a.out, index=False)
        print("written:", a.out, file=sys.stderr)


if __name__ == "__main__":
    main()
