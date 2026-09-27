# Tested software versions & pinned parameters

This file records the **exact toolchain the EpiGuide / NeuroScore workflow was developed
and tested with**, so the pipeline is reproducible. Third-party command-line tools
(especially `modkit`) have changed their flags and output between releases — please match
the versions below, or read the compatibility notes.

## Core runtime (Python, pip)

| Package | Tested version | Constraint in `requirements.txt` |
|---|---|---|
| Python | 3.12 | 3.10+ |
| numpy | 1.26 | `numpy<2` (required by torch 2.2) |
| torch | 2.2 | `torch>=2.0,<2.3` |
| torch-geometric | 2.8 | `torch-geometric>=2.3` |
| pandas | 2.1 | `pandas>=1.5` |

## External tools (not installed by pip)

| Tool | **Tested / pinned version** | Notes |
|---|---|---|
| **modkit** | **0.3.0** | Used for `bam_to_calls.py`. See command + compatibility below. |
| samtools | 1.17 | sort/index/fastq for alignment (upstream of this repo) |
| minimap2 | 2.28-r1221 | alignment (upstream of this repo) |

Install the tested `modkit` (e.g. via bioconda):

```bash
conda install -c bioconda ont-modkit=0.3.0
# or download the 0.3.0 release binary:
# https://github.com/nanoporetech/modkit/releases/tag/v0.3.0
```

Verify:

```bash
modkit --version   # expect: mod_kit 0.3.0
```

## `modkit pileup` — exact tested command

`bam_to_calls.py` runs, per sample:

```bash
modkit pileup <input.bam> <pileup.bed> \
    --include-bed model/model_cpgs_hg38.bed \
    --threads 4 \
    [--ref hg38.fa]        # optional
```

The BAM must contain ONT modified-base tags (`MM`/`ML`, 5mC in CpG context) and be aligned
to **GRCh38 / hg38**. Output rows with `mod_code == "m"` (5mC) are kept and, per CpG,
`methylation_calls = Nmod`, `unmethylation_calls = Ncanonical`,
`total_calls = Nmod + Ncanonical`.

## Compatibility with newer modkit (≥ 0.5, incl. 0.6.3)

In modkit **0.3.0** the command above runs as-is. In **newer releases the `--include-bed`
region filter can require an explicit motif/base specification** (e.g. it may exit with
`--include-bed requires either --motif or --modified-bases`), which caused reproduction
failures reported by reviewers.

Because `bam_to_calls.py` **already restricts the pileup to the panel CpGs in software**
(the pileup is merged against `model/model_cpgs_hg38.bed` by genomic coordinate), the
`--include-bed` flag is only a speed optimisation and is **not required for correctness**.
The code therefore:

1. tries the region-restricted command above; and
2. if `modkit` exits non-zero, **automatically retries an unrestricted**
   `modkit pileup <bam> <out> --threads N [--ref ...]` and filters to the panel in pandas.

The retry produces **identical NeuroScore results**; it is only slower (genome-wide pileup).

If you prefer to keep region+motif restriction on a newer modkit, the equivalent modern
invocation is:

```bash
modkit pileup <input.bam> <pileup.bed> \
    --include-bed model/model_cpgs_hg38.bed \
    --cpg --ref hg38.fa \
    --threads 4
```

(`--cpg` needs `--ref`.) The downstream parsing in `bam_to_calls.py` is unchanged.

## Provenance of the analysed cohort (for reference)

The study BAMs were basecalled with `dna_r10.4.1_e8.2_400bps_hac@v4.3.0` (R10.4.1, HAC),
aligned with minimap2 2.28-r1221 (`-x map-ont -a -y --secondary=no`) to `hg38.fa.gz`, and
piled up with modkit 0.3.0. crossNN/nanoDx classification ran through ROBIN 0.5.
