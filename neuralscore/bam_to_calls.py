"""BAM -> per-CpG methylation calls at the model's CpG sites, via modkit.

Requires `modkit` (https://github.com/nanoporetech/modkit) on PATH and a BAM with
MM/ML modified-base tags (ONT 5mC), aligned to GRCh38/hg38. Produces a table with
ID_REF, methylation_calls, unmethylation_calls, total_calls for predict_from_calls().
"""
import os, shutil, subprocess, tempfile
import pandas as pd

# ONT bedMethyl columns (modkit pileup)
_BM_COLS = ["chrom","start","end","mod_code","score","strand","tstart","tend","color",
            "Nvalid","fraction","Nmod","Ncanonical","Nother","Ndelete","Nfail","Ndiff","Nnocall"]


def bam_to_calls(bam, model_bed, out_tsv=None, modkit="modkit", threads=4,
                 mod_code="m", ref=None):
    if shutil.which(modkit) is None:
        raise RuntimeError(
            f"'{modkit}' not found on PATH. Install modkit "
            "(https://github.com/nanoporetech/modkit) or pass --calls with a precomputed table.")

    bed = pd.read_csv(model_bed, sep="\t", header=None, names=["chrom","start","end","ID_REF"])
    bed["chrom"] = bed["chrom"].astype(str)

    with tempfile.TemporaryDirectory() as tmp:
        pileup = os.path.join(tmp, "pileup.bed")
        cmd = [modkit, "pileup", bam, pileup, "--include-bed", model_bed,
               "--threads", str(threads)]
        if ref:
            cmd += ["--ref", ref]
        subprocess.run(cmd, check=True)
        bm = pd.read_csv(pileup, sep=r"\s+", header=None, names=_BM_COLS, engine="python")

    bm = bm[bm["mod_code"] == mod_code].copy()
    bm["chrom"] = bm["chrom"].astype(str)
    merged = bm.merge(bed, on=["chrom", "start"], how="inner")
    calls = pd.DataFrame({
        "ID_REF": merged["ID_REF"],
        "methylation_calls": merged["Nmod"].astype(int),
        "unmethylation_calls": merged["Ncanonical"].astype(int),
        "total_calls": (merged["Nmod"] + merged["Ncanonical"]).astype(int),
    })
    calls = calls.groupby("ID_REF", as_index=False).sum()
    if out_tsv:
        calls.to_csv(out_tsv, sep="\t", index=False)
    return calls
