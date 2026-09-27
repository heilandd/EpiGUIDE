"""BAM -> per-CpG methylation calls at the model's CpG sites, via modkit.

Requires `modkit` (https://github.com/nanoporetech/modkit) on PATH and a BAM with
MM/ML modified-base tags (ONT 5mC), aligned to GRCh38/hg38. Produces a table with
ID_REF, methylation_calls, unmethylation_calls, total_calls for predict_from_calls().

Tested with modkit 0.3.0 (see VERSIONS.md). modkit's command-line interface has changed
between releases: in newer versions (>= 0.5, incl. 0.6.3) the `--include-bed` region
filter can require an explicit `--motif`/`--modified-bases` specification and otherwise
exits non-zero. Because this module already restricts the pileup to the panel CpGs in
software (merge against `model_cpgs_hg38.bed`), `--include-bed` is only a speed
optimisation; if the region-restricted call fails, we automatically retry an unrestricted
`modkit pileup` and filter in pandas, giving identical results.
"""
import os, re, shutil, subprocess, tempfile, warnings
import pandas as pd

# ONT bedMethyl columns (modkit pileup)
_BM_COLS = ["chrom","start","end","mod_code","score","strand","tstart","tend","color",
            "Nvalid","fraction","Nmod","Ncanonical","Nother","Ndelete","Nfail","Ndiff","Nnocall"]

TESTED_MODKIT = "0.3.0"


def modkit_version(modkit="modkit"):
    """Return modkit's version string (e.g. '0.3.0'), or None if it can't be determined."""
    try:
        out = subprocess.run([modkit, "--version"], capture_output=True, text=True, check=True)
        m = re.search(r"(\d+\.\d+\.\d+)", (out.stdout or "") + (out.stderr or ""))
        return m.group(1) if m else None
    except Exception:
        return None


def bam_to_calls(bam, model_bed, out_tsv=None, modkit="modkit", threads=4,
                 mod_code="m", ref=None):
    if shutil.which(modkit) is None:
        raise RuntimeError(
            f"'{modkit}' not found on PATH. Install the tested version "
            f"(modkit {TESTED_MODKIT}; see VERSIONS.md) or pass --calls with a "
            "precomputed table.")

    ver = modkit_version(modkit)
    if ver and ver != TESTED_MODKIT:
        warnings.warn(
            f"modkit {ver} detected; this workflow was tested with modkit {TESTED_MODKIT} "
            "(see VERSIONS.md). Command-line options differ between releases; the "
            "region-restricted call will be retried unrestricted if it fails.")

    bed = pd.read_csv(model_bed, sep="\t", header=None, names=["chrom","start","end","ID_REF"])
    bed["chrom"] = bed["chrom"].astype(str)

    with tempfile.TemporaryDirectory() as tmp:
        pileup = os.path.join(tmp, "pileup.bed")
        base = [modkit, "pileup", bam, pileup, "--threads", str(threads)]
        if ref:
            base += ["--ref", ref]
        # 1) tested command: restrict pileup to panel regions (fast).
        restricted = base + ["--include-bed", model_bed]
        try:
            subprocess.run(restricted, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            # 2) fallback for modkit versions where --include-bed needs a motif/base spec:
            #    run genome-wide pileup; the panel filter is applied below in software.
            warnings.warn(
                "modkit pileup with --include-bed failed (likely a modkit-version CLI "
                f"change; see VERSIONS.md). stderr:\n{(e.stderr or '').strip()}\n"
                "Retrying an unrestricted 'modkit pileup' and filtering to the panel in "
                "software (identical result, slower).")
            subprocess.run(base, check=True)
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
