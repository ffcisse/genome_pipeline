<a id="readme-top"></a>

<div align="center">

# Genome Comparison Pipeline

**A reproducible Snakemake pipeline for comparative genomics — protein/CDS property computation, cross-genome statistics, static figures, and an interactive dashboard, all driven by two config files.**

[![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A59.23.1-039475?style=for-the-badge&logo=snakemake&logoColor=white)](https://snakemake.readthedocs.io/)
[![Python](https://img.shields.io/badge/python-3.11-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org/)

</div>

<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#about-the-project">About The Project</a>
      <ul><li><a href="#built-with">Built With</a></li></ul>
    </li>
    <li><a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a>
      <ul>
        <li><a href="#configuration">Configuration</a></li>
        <li><a href="#input-requirements">Input Requirements</a></li>
        <li><a href="#running-it">Running It</a></li>
        <li><a href="#dependencies--environments">Dependencies &amp; Environments</a></li>
        <li><a href="#outputs">Outputs</a></li>
        <li><a href="#visualizations-phase-5">Visualizations (Phase 5)</a></li>
        <li><a href="#interactive-dashboard-phase-6a">Interactive Dashboard (Phase 6a)</a></li>
        <li><a href="#interpreting-the-statistics">Interpreting the Statistics</a></li>
      </ul>
    </li>
    <li><a href="#known-limitations">Known Limitations</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>

## About The Project

Point this pipeline at a set of genomes' protein + CDS FASTA files and it computes per-protein and
per-gene properties, merges them into cross-genome master tables, runs configurable
group-comparison statistics (effect sizes and a leave-one-subgroup-out sensitivity analysis), and
produces both a full set of static figures and a single-file interactive dashboard for exploring
the results.

**In:** one protein FASTA + one CDS FASTA per genome, plus a genome metadata table
(`config/genomes.tsv`) with at least a genome ID and two grouping columns.

**Out:**
- Per-genome tables: physicochemical protein properties, intrinsic disorder, codon usage/ENC/GC
  content
- Cross-genome master tables (all genomes' per-gene data joined with group labels)
- A species-level summary table (median/mean/std of every property, per genome)
- Effect-size tables comparing genomes by any grouping column you configure (e.g. lifestyle,
  lineage, or your own)
- A leave-one-subgroup-out sensitivity analysis, to check whether an apparent group difference is
  actually being driven by one phylogenetic subgroup rather than the grouping you think you're
  testing
- 457 static figures (boxplots, violins, distributions, PCA, clustering, effect-size forest plots,
  a sensitivity heatmap), as both PNG and PDF
- A single standalone HTML dashboard for interactively exploring the same properties/groupings

### Built With

[![Snakemake](https://img.shields.io/badge/Snakemake-039475?style=for-the-badge&logo=snakemake&logoColor=white)](https://snakemake.readthedocs.io/)
[![Python](https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org/)
[![pandas](https://img.shields.io/badge/pandas-150458?style=for-the-badge&logo=pandas&logoColor=white)](https://pandas.pydata.org/)
[![Biopython](https://img.shields.io/badge/Biopython-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://biopython.org/)
[![metapredict](https://img.shields.io/badge/metapredict-orange?style=for-the-badge)](https://github.com/idptools/metapredict)
[![codonW](https://img.shields.io/badge/codonW-1.4.4-lightgrey?style=for-the-badge)](https://anaconda.org/bioconda/codonw)
[![Plotly.js](https://img.shields.io/badge/Plotly.js-3F4F75?style=for-the-badge&logo=plotly&logoColor=white)](https://plotly.com/javascript/)
[![SLURM](https://img.shields.io/badge/SLURM-blue?style=for-the-badge)](https://slurm.schedmd.com/)

| Tool | Used for |
|---|---|
| **Snakemake** | Pipeline orchestration — every rule, wildcard, and DAG dependency in `workflow/` |
| **Python 3.11** | Every rule's implementation (`workflow/scripts/*.py`) |
| **pandas / numpy / scipy** | Tabular data, statistics (Mann-Whitney U, effect sizes) |
| **Biopython** | FASTA parsing (`qc`, `parse`) |
| **metapredict** | Intrinsic disorder prediction (Phase 2b, PyTorch-based) |
| **codonW** | ENC, GC, GC3s, codon usage (Phase 3, standalone C program) |
| **Plotly.js** | The interactive dashboard's plots (Phase 6a, vendored/inlined, not a CDN) |
| **SLURM** | Optional cluster execution on DORI — see [Running It](#running-it) |

<p align="right">(<a href="#readme-top">back to top</a>)</p>


## Project Origin

This pipeline was developed and validated using a comparative genomics project on **Rhodophyta
(red algae)**, comparing extremophilic and mesophilic lineages. While the pipeline itself is
genome-agnostic (see [Configuration](#configuration)), the dataset below is the one it was built
and tested against, and is used throughout this README for examples.

**Study design:** 9 genomes across two lifestyle groups and three phylogenetic lineages.

| genome | lifestyle | lineage |
|---|---|---|
| Cyamer1 | extremophile | Cyanidiales |
| CyamerSoos_1_1 | extremophile | Cyanidiales |
| Cyanyang1 | extremophile | Cyanidiales |
| Galph1_1 | extremophile | Galdieria |
| Galsul1 | extremophile | Galdieria |
| Galyel1 | extremophile | Galdieria |
| Porcrue1 | mesophile | Mesophile_lineage |
| Porpu1328_1 | mesophile | Mesophile_lineage |
| Rhomari1 | mesophile | Mesophile_lineage |


<p align="right">(<a href="#readme-top">back to top</a>)</p>


## Getting Started

### Prerequisites

- **conda or mamba**, to create the Snakemake environment and every per-rule environment
  (`--use-conda` activates them automatically — nothing to install by hand beyond conda itself).
- Enough local cores for a quick test, or access to a **SLURM cluster** for a real run (SLURM is
  optional — see [Running It](#running-it) for both paths).
- Your genomes' protein + CDS FASTA files (or a way to fetch them — see
  [Input Requirements](#input-requirements)).

### Installation

```bash
git clone git@github.com:ffcisse/genome_pipeline.git
cd genome_pipeline

# 1. Conda env with Snakemake itself
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake

# 2. Point config at your genomes
#    Open each file below in a text editor (e.g. `vim`, `nano`, or `code`) and fill in your
#    own genome names, paths, and grouping labels.
#
#    - config/genomes.tsv: one row per genome (genome ID, input paths, and any grouping
#      columns you want to use for comparisons -- e.g. lifestyle, lineage)
#    - config/config.yaml: set input.protein_dir / input.cds_dir (or staging.source_dir),
#      plus which genomes.tsv columns to use as your primary grouping and subgroup
#
#    See [Configuration](#configuration) for the full column reference and a worked example.
vim config/genomes.tsv
vim config/config.yaml
```

No code editing required to point this at a different genome set — just `config/genomes.tsv` and
`config/config.yaml`, including the grouping columns used for effect sizes/sensitivity analysis
(see [Configuration](#configuration)). See [Known Limitations](#known-limitations).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Usage

```bash
# Dry run -- see what would happen, without running anything
snakemake -n

# Run it
snakemake --cores 4 --use-conda          # locally (small datasets / a login node's own cores)
# or, on a SLURM cluster:
workflow/scripts/submit_phase1.sh        # see Running It below for the full phase list
```

`rule all`'s default target is the dashboard itself (`results/dashboard/proteome_dashboard.html`),
so a plain `snakemake --cores 4 --use-conda` (or the last submit script in the phase list, since
each one just runs bare `snakemake`) pulls in every phase through Phase 6B, figures and dashboard
included — see [Interactive Dashboard](#interactive-dashboard-phase-6a) for what to do with the
result.

### Configuration

Edit these two files to run on a new genome set.

#### `config/config.yaml`

| Key | Required? | Meaning |
|---|---|---|
| `input.mode` | yes | Only `existing_dir` is implemented. |
| `input.protein_dir` | yes | Directory of per-genome protein FASTA, named `<genome_id>.fasta` (or see `staging` below). |
| `input.cds_dir` | yes | Same, for CDS FASTA. |
| `staging.source_dir` | no | Set this instead of pre-flattening `protein_dir`/`cds_dir` yourself — see [Input Requirements](#input-requirements). Omit the whole `staging:` block if your files are already flat. |
| `staging.protein_subdir` | no (default `proteome_files`) | Subdirectory name under `source_dir` holding protein downloads. |
| `staging.cds_subdir` | no (default `cds_files`) | Same, for CDS downloads. |
| `genome_table` | yes | Path to the genome metadata TSV (see below). Default `config/genomes.tsv`. |
| `output_dir` | yes | Where results land. Default `results`. |
| `sensitivity.primary_grouping` | yes | Column in `genome_table` to test as the "main" grouping (e.g. `lifestyle`) — must have exactly 2 distinct values. |
| `sensitivity.subgroup_column` | yes | Column in `genome_table` whose values get dropped one at a time in the sensitivity sweep (e.g. `lineage`). |
| `group_value_order` | no | Optional: pins which value of a grouping column is "A" (vs "B") in `cles`/`rank_biserial`, and the order pairwise comparisons are generated in for a >2-value column. A dict of `column: [value, value, ...]`; omit a column (or the whole block) for alphabetical order instead. |
| `dashboard.sample_per_genome` / `dashboard.sample_seed` | no (defaults `2500` / `42`) | How many proteins per genome the dashboard's distribution views sample, and the fixed seed for reproducibility. See [Interactive Dashboard](#interactive-dashboard-phase-6a). |
| `slurm.account` / `slurm.qos` / `slurm.partition` / `slurm.mail_user` | no, but **deployment-specific** | Passed to `sbatch` by the `submit_phase*.sh` wrappers. **You must change these** — they're this deployment's cluster allocation and email, not yours. Leave `partition: ""` if your cluster doesn't need one. |

#### `config/genomes.tsv`

Tab-separated, one row per genome:

| Column | Required? | Meaning |
|---|---|---|
| `genome_id` | **yes** | Must match the `<genome_id>` in your FASTA filenames (or staging source tree). Drives every `{genome}` wildcard in the pipeline. |
| `name` | no | Free-text display name, not used in any computation. |
| *(two grouping columns)* | **yes** | Column names and values are fully configurable — see below. |

**Groupings come from `config/genomes.tsv`, not from code.** Give each genome's row whatever label
columns describe your study, then tell the pipeline which two of those columns to use for the
built-in statistics via `config.yaml`'s `sensitivity.primary_grouping`/`subgroup_column`:

- **Primary grouping** — your main, 2-group comparison (e.g. `lifestyle`: extremophile vs.
  mesophile). This is the "A vs. B" every effect-size table tests.
- **Subgroup column** — a finer-grained column nested inside the primary grouping (e.g.
  `lineage`, a phylogenetic clade). The leave-one-subgroup-out sensitivity analysis drops one
  subgroup value at a time, to check whether the primary-grouping effect is really just one
  subgroup driving it rather than the thing you think you're measuring.

A tiny example `config/genomes.tsv`:

| genome_id | lifestyle | lineage |
|---|---|---|
| Genome_A | extremophile | Clade_1 |
| Genome_B | extremophile | Clade_1 |
| Genome_C | mesophile | Clade_2 |

Column names and values here are entirely up to you — nothing in the code is hardcoded to
`lifestyle`/`lineage`.

**Which value is "A"?** For any grouping column, `cles`/`rank_biserial` need a direction — which
value counts as "A" (vs "B"), and, for a >2-value column, what order pairwise comparisons are
generated in. That's an analyst framing choice (which group are you testing *for*), not something
derivable from the data alone, so it's a separate, optional config key:
`config.yaml`'s `group_value_order`. Leave a column out of it (or omit the whole block) and its
values sort alphabetically instead — either way, nothing needs a code change.

**A third grouping dimension** (beyond the two `sensitivity:` already names) isn't wired up by
default — `workflow/Snakefile`'s `GROUPING_COLUMNS` list is built from exactly
`sensitivity.primary_grouping`/`subgroup_column`, so those are the only two columns that
automatically flow into the master tables and get an `effect_sizes_<grouping>.csv`. Every script
already accepts an arbitrary set of group columns (nothing to edit there); adding a third means
adding a config key and one line in the Snakefile to fold it into `GROUPING_COLUMNS` — see
[Known Limitations](#known-limitations) for exactly what that would look like.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Input Requirements

Each genome needs a protein FASTA and a CDS FASTA. Two ways to point the pipeline at them:

1. **Flat layout (simplest):** put `<genome_id>.fasta` (optionally `.gz` or `.tar.gz` — the
   `qc`/`parse` scripts handle plain, gzipped, and gzipped-tar-with-one-member transparently) in
   `input.protein_dir` / `input.cds_dir`.
2. **Nested download tree:** if your files live several directories deep (e.g. a Mycocosm-style
   download), set `staging.source_dir` in `config.yaml`. The `stage_inputs` rule finds, for each
   genome, the single file matching `source_dir/<protein_subdir|cds_subdir>/<genome_id>/**/*.fasta*`
   and symlinks it into the flat layout above — real extension preserved, no decompression.

**Optional JGI download helper:** `workflow/scripts/download_from_jgi.sh` documents how to fetch
FASTA from JGI Mycocosm/Phycocosm. Be aware:
- It is **not wired into the pipeline** — a standalone, manual step.
- It must run on a **NERSC Perlmutter DTN** (JGI's download service isn't reachable from most
  other places), not wherever you run the rest of this pipeline.
- **As shipped, it's a documented stub** — it prints a message and exits immediately; the real
  `portal-apps-bootstrap.sh` invocation is commented out below the stub, to be filled in once you
  have that tool available. Read the comments in the script before relying on it.
- If your genomes are on **private** Mycocosm/Phycocosm portals, the real invocation needs
  `--use-non-public-portals true` — omitting it silently returns nothing for private-portal
  genomes.
- None of this matters if you already have the FASTA files from anywhere else — just point
  `config.yaml` at them.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Running It

#### Dry run

```bash
snakemake -n
```

Shows the full job plan without executing anything: which rules will run, for which genomes, and
why (missing output, changed input, etc.). Always run this after a config change. Add `-p` to
also print each job's actual shell command, or `--dag | dot -Tpng > dag.png` to render the
dependency graph visually.

#### Per-phase SLURM submission

Each phase has a `submit_phaseN.sh` wrapper (reads `config.yaml`'s `slurm:` block and passes it to
`sbatch`) around a `run_phaseN.sbatch` script:

| Script | Runs | Resources | Notes |
|---|---|---|---|
| `submit_phase1.sh` | `stage_inputs`→`qc`→`parse` | 4 CPU / 16G | Light |
| `submit_phase2a.sh` | `protein_properties` | 4 CPU / 16G | Light |
| `submit_phase2b.sh` | `disorder` | **64 CPU (exclusive node) | **Heavy — see warning below** |
| `submit_phase3.sh` | `cds_properties` | 4 CPU / 16G | Light |
| `submit_phase4.sh` | `merge_*`/`species_summary`/`effect_sizes`/`sensitivity_*` | 4 CPU / 16G | Light |
| `submit_phase5.sh` | `plot_*` (457 figures) | 4 CPU / 16G | Light |
| `submit_phase6.sh` | `interactive_dashboard_*` | 4 CPU / 16G | Light |

Each script runs `snakemake --cores $SLURM_CPUS_PER_TASK --use-conda` locally inside one SLURM
allocation (not Snakemake's separate cluster-executor-plugin model), so running an earlier
phase's script again after later work is done is harmless — it just confirms everything's already
up to date. Phases 1, 2a, 3, 4, and 5 build their outputs as part of `rule all`'s default target,
so their submit scripts run bare `snakemake`; Phase 2b's `disorder` output is deliberately *not*
part of `rule all` (see below), so `submit_phase2b.sh` names its targets explicitly.

**How the phase scripts actually work.** Snakemake resolves dependencies backward from a target,
not forward from a starting phase. Each `submit_phaseN.sh` script just tells Snakemake "build this
output file" — Snakemake then automatically figures out everything upstream that's needed to
produce it, and runs *only* whatever is missing or out of date.

This means:
- **You don't have to run phases in order.** If you run `submit_phase6.sh` on a fresh checkout with
  nothing done yet, it will run the entire pipeline end-to-end (staging through the dashboard) in
  one submission — there's no requirement to run phases 1-5 first.
- **Running a later phase after earlier ones are done only builds what's new.** If phases 1-5 are
  already complete, `submit_phase6.sh` will skip everything that already exists and just build the
  summaries/visuals/dashboard steps still missing.
- **The phase-numbered scripts exist for convenience and incremental verification** (so you can
  check each stage's output before moving on) Feel free to jump straight to whichever `submit_phaseN.sh` produces the output you
  actually want.

Submit via the wrapper, not `sbatch run_phaseN.sbatch` directly, unless your site has defaults for
account/QOS/mail-user — the wrapper is what supplies those from `config.yaml`.

> [!WARNING]
> **Phase 2b (`disorder`) is heavy.** It loads a PyTorch model and runs inference over every
> protein. Manual single-threaded testing needed ~2h/genome, which is why this rule requests a
> full 64-core exclusive node (the resource sizing behind `run_phase2b.sbatch`'s 1.5h budget); the
> pipeline's actual batched implementation is considerably faster in practice — the verified
> 9-genome, 61,349-protein run finished in under 10 minutes on that allocation. Budget for the
> heavier estimate regardless (node contention and dataset size both vary), and don't run this
> rule inline on a shared login node.

**Logs** land in `logs/phaseN_<jobid>.{out,err}`.

**SLURM values you must change:** `config.yaml`'s `slurm:` block (`account`, `qos`, `partition`,
`mail_user`) is this deployment's cluster allocation and email address — update it for yours
before submitting anything. See [`SLURM.md`](SLURM.md) for DORI/SLURM-specific notes beyond what's
in this README (the `--use-conda`-in-SLURM testing discipline, per-rule `resources:`, and why
there's no cluster-executor-plugin model here).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Dependencies & Environments

Every rule that needs real dependencies declares its own conda env under `workflow/envs/`,
activated automatically via `--use-conda`:

| Env | Used by | Key packages |
|---|---|---|
| `analysis.yaml` | `qc`, `parse`, `protein_properties` | biopython, pandas, numpy, scipy |
| `disorder.yaml` | `disorder` | metapredict (pip; pulls PyTorch) |
| `codonw.yaml` | `cds_properties` | codonW 1.4.4 (bioconda) |
| `summaries.yaml` | Phase 4 rules, Phase 6a's `dashboard_data`/`dashboard` | pandas, numpy, scipy |
| `visuals.yaml` | Phase 5's `plot_*` rules | pandas, numpy, scipy, scikit-learn (PCA), matplotlib, seaborn |

- **codonW** comes from bioconda (`codonw=1.4.4=h7b50bb2_7`, pinned to a build already validated
  on this deployment) — it's a standalone C program, not a Python package.
- **metapredict** comes from PyPI via the env's `pip:` section (not on conda-forge/bioconda) and
  pulls in PyTorch as a dependency — this is why `disorder.yaml` is CPU-inference-heavy to
  install(see the SLURM warning above).

> [!WARNING]
> **Known issue:** `--use-conda` environments have broken *inside SLURM
> jobs* on this deployment more than once — while working fine when tested interactively on a
> login node. Root causes varied (a cold/first-touch environment on a freshly allocated node in
> one case; a genuine `libstdc++` ABI conflict between pip-installed PyTorch and conda-forge-built
> numpy/scipy in another — see `workflow/rules/disorder.smk`'s shell-command comment for the full
> story and fix). **Test any new or modified conda env inside a real SLURM job before trusting
> it** — a login-node pass is not sufficient evidence it'll work under `sbatch`.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Outputs

```
results/
├── qc/<genome>.qc.done                          # human-readable QC report
├── parsed/<genome>/{protein,cds}_table.csv       # canonical per-gene sequence tables
├── protein_properties/<genome>/protein_properties.csv
├── disorder/<genome>/disorder.csv
├── cds_properties/<genome>/{cds_properties,codon_usage}.csv
├── summaries/
│   ├── master_protein_table.csv                 # all genomes' protein_properties + disorder + group labels
│   ├── master_cds_table.csv                     # all genomes' cds_properties + codon_usage + group labels
│   ├── species_summary.csv                      # one row per genome
│   ├── effect_sizes_<grouping>.csv               # one per configured grouping (default: lifestyle, lineage)
│   ├── sensitivity_drop_<subgroup>.csv           # one per subgroup value
│   └── sensitivity_leave_one_out.csv             # all subgroups combined
├── plots/                                        # Phase 5 -- 457 static figures, PNG (300dpi) + PDF each
│   ├── boxplots/            # box+violin per property, by primary grouping / subgroup / per-species
│   ├── distributions/       # histogram+KDE per property, by primary grouping / subgroup
│   ├── cds_distributions/   # same as distributions/, for CDS-level properties (ENC/GC/GC3s/length)
│   ├── pca/                 # per-protein and per-species PCA, colored by each grouping, + PC1/PC2 loadings
│   ├── clustering/          # genome x property hierarchical clustermap; property-property Spearman clustermap
│   ├── effect_sizes/        # forest plots of rank-biserial, sorted by |magnitude|
│   └── sensitivity/         # leave-one-out shrinkage heatmap
└── dashboard/                                    # Phase 6a
    ├── data.json                # the JSON payload embedded in the dashboard (also useful standalone)
    └── proteome_dashboard.html  # the standalone interactive dashboard -- see below
```

<details>
<summary><strong>Key table columns</strong> (click to expand)</summary>

**`master_protein_table.csv`** (29 cols; one row per protein, 61,349 in the verified dataset) —
`genome`, `protein_id`, `sequence`, `length`, composition (`pct_charged`, `pct_acidic`,
`pct_basic`, `pct_hydrophobic`, `pct_polar`, `pct_aromatic`, `pct_special`), `pI`, `gravy`,
`instability_index`, `net_charge_pH7`, `aliphatic_index`, `thermostable_fraction`,
`cysteine_fraction`, `carbon_oxidation_state`, `charge_density`, aggregation
(`agg_mean_a3v`, `agg_Na4vSS`, `agg_hotspot_fraction`), disorder (`disorder_mean`,
`disorder_fraction`, `longest_idr`, `n_idrs`), plus your configured group columns
(`lifestyle`, `lineage` by default).

**`master_cds_table.csv`** (77 cols) — `genome`, `cds_id`, `length`, `is_multiple_of_3`, `enc`
(effective number of codons, ~20-61), `gc`, `gc3s` (both proportions in **[0, 1]**, not
percentages), `start_codon`, `stop_codon`, `has_canonical_start`, `has_canonical_stop`, 64
`codon_<TRIPLET>` raw-count columns, plus group columns.

**`species_summary.csv`** (one row per genome, 290 cols in the verified dataset) — `genome`,
group columns, `n_proteins`, `n_cds`, then `<property>_median`/`_mean`/`_std` for every numeric
property in both master tables (95 properties × 3 stats).

**`effect_sizes_<grouping>.csv`** — `property`, `table` (`protein`/`cds` — needed because a few
property names, e.g. `length`, exist in both tables with different meanings), `group_a`,
`group_b`, `n_a`, `n_b`, `median_a`, `median_b`, `p_value`, `cles`, `rank_biserial`. One row per
property for a 2-value grouping; one row per property per pair for an N-value grouping (e.g. 3
pairwise rows per property for a 3-value `lineage`). Sorted by `|rank_biserial|` descending.

**`sensitivity_leave_one_out.csv`** — `excluded_subgroup`, `property`, `table`,
`rank_biserial_full`, `rank_biserial_excluded`, `shrinkage`. See
[Interpreting the Statistics](#interpreting-the-statistics) for what `shrinkage` means.

</details>

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Visualizations (Phase 5)

`results/plots/` holds 457 static figures (PNG at 300dpi + PDF, every figure in both formats),
generated by the `plot_*.py` scripts in `workflow/rules/visuals.smk` directly from Phase 4's
already-computed tables — no statistics are recomputed for plotting.

| Subdirectory | Contents |
|---|---|
| `boxplots/` | Box+violin plots of every protein property, by primary grouping, by subgroup, and by species (median-labeled, colored by subgroup) |
| `distributions/` | Histogram + KDE of every protein property, by primary grouping and by subgroup |
| `cds_distributions/` | Same as `distributions/`, for CDS-level scalar properties (ENC, GC, GC3s, length, start/stop codon) |
| `pca/` | PCA of the protein property matrix at per-protein and per-species (genome-median) resolution, colored by each grouping, plus PC1/PC2 loadings |
| `clustering/` | Hierarchical clustermap of genome × property (z-scored medians), and a property-property Spearman correlation clustermap |
| `effect_sizes/` | Forest/bar plots of rank-biserial effect sizes, sorted by `\|magnitude\|` — one for the primary grouping, one per pairwise subgroup comparison |
| `sensitivity/` | The leave-one-out: a heatmap of shrinkage (top properties × excluded subgroup) |

All of it excludes `cds_properties.py`'s 64 raw `codon_<TRIPLET>` count columns from per-property
figures and effect-size/sensitivity rankings (they'd otherwise dominate a combined ranking purely
by count) — `config.yaml`'s `visuals.exclude_properties` can drop additional specific properties
from the per-property figures, and `visuals.top_n_effect_sizes` caps how many properties the
effect-size/sensitivity figures show (full results always remain in the Phase 4 CSVs regardless).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Interactive Dashboard (Phase 6)

`results/dashboard/proteome_dashboard.html` is **one standalone HTML file** (~5MB for the real
9-genome dataset) that opens directly in a browser: Plotly.js is embedded rather than loaded from a CDN, and the
data it plots is embedded JSON. Built by two rules (`dashboard_data` → `dashboard`;
see `workflow/rules/dashboard.smk`), part of `rule all`'s default target.

**Getting it off the cluster:** it's a single file, so a plain `scp` works —

```bash
scp <your-username>@<cluster-host>:/path/to/genome_pipeline/results/dashboard/proteome_dashboard.html ~/Desktop/
```

then double-click it (or open it from your browser's File → Open) on your laptop. No further setup.

**What's in it** — five sections, navigated from the left sidebar:
- **Overview** — genome/group counts, total proteins/CDS, and a Phase 1 QC summary table.
- **Property Explorer** — pick a property, a grouping (species / primary grouping / subgroup), and
  a plot type (box, violin, histogram, or density/KDE); the plot updates live. This is the
  interactive replacement for paging through Phase 5's per-property static figures.
- **Species View** — a sortable per-genome summary table (click a column header to sort), plus the
  same property/plot-type view scoped to genomes.
- **Effect Sizes** — an interactive forest plot of rank-biserial effect sizes, switchable between
  the primary grouping and any subgroup pairwise comparison, sorted by magnitude. Clicking a bar
  jumps to that property in Property Explorer, including CDS-level properties (Phase 6b).
- **Sensitivity** — the leave-one-subgroup-out shrinkage heatmap, with a plain-language explanation
  of positive vs. negative shrinkage (same interpretation as
  [Interpreting the Statistics](#interpreting-the-statistics)).

> [!NOTE]
> **The sampling caveat** (read this before trusting a shape you see in the dashboard): every
> summary number in the dashboard — medians, the *exact* box-plot quartiles, effect sizes,
> sensitivity — is computed from the **full** dataset
> , never recomputed in the browser. But the violin/histogram/density views, and the
> Species View comparison plot when set to those modes, are drawn from a **downsampled**
> per-protein sample: up to 2,500 proteins per genome, fixed seed (42) for reproducibility,
> ~22,500 rows total for the real dataset. The UI captions every sampled view accordingly
> ("distributions shown from a sample..."). If you need an exact distribution shape rather than a
> representative one, use the box plot view or the underlying CSVs directly.


<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Interpreting the Statistics

**p-values don't indicate much.** With ~61,000 proteins/CDS in the example, Mann-Whitney p-values collapse
to ~0 for almost every property regardless of whether the difference is biologically meaningful —
they're reported for reference, not as evidence of importance.

**Use CLES and rank-biserial instead:**
- `cles` (common-language effect size) = the probability a randomly picked member of `group_a`
  exceeds a randomly picked member of `group_b`. 0.5 = no difference.
- `rank_biserial` = `2*cles - 1`, rescaled to **-1..+1** (0 = no difference). Rough magnitude
  guide: `|0.1|` small, `|0.3|` medium, `|0.5|` large.

**Why the leave-one-subgroup-out sensitivity analysis exists:** when your subgroup column is
nested inside your primary grouping (e.g. every genome's `lineage` implies its `lifestyle`), an
apparent primary-grouping effect can actually be driven entirely by one subgroup, not by the
grouping variable itself — a phylogenetic confound rather than a real effect of what you think
you're measuring. `sensitivity_leave_one_out.csv` quantifies each subgroup's individual
contribution by dropping it and recomputing:
- **`shrinkage` > 0**: the excluded subgroup was inflating the apparent effect — remove it and the
  primary grouping separates less. The bigger the shrinkage, the more that subgroup alone was
  driving the signal.
- **`shrinkage` < 0**: the opposite — removing that subgroup made the remaining groups separate
  *more*, meaning it had been diluting the effect (likely because it overlaps the other primary
  group rather than being distinct from it).
- **`shrinkage` ≈ 0**: the property's effect doesn't depend much on that particular subgroup.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Known Limitations

- **Only two grouping columns are wired up by default.** `workflow/Snakefile`'s `GROUPING_COLUMNS`
  is built from exactly `sensitivity.primary_grouping`/`subgroup_column` — those are the only two
  columns that automatically get merged into the master tables and get their own
  `effect_sizes_<grouping>.csv`. Adding a third grouping dimension doesn't need any script changes
  (every script already accepts an arbitrary set of group columns via `--group-columns`/
  `--genomes-tsv`) — just a new config.yaml key and one line in the Snakefile appending it to
  `GROUPING_COLUMNS`.
- **`download_from_jgi.sh` is a documented stub**, not a working downloader (see
  [Input Requirements](#input-requirements)) — it requires a JGI-provided tool
  (`portal-apps-bootstrap.sh`) this repo doesn't include. Prints a message and exits immediately;
  read the comments in the script for the real invocation, commented out below the stub.
- **SLURM account/QOS/partition/email are this deployment's values** and must be changed in
  `config.yaml` by anyone else running this pipeline.
- **`staging.source_dir`** (in the shipped `config.yaml`) is an absolute path on this deployment's
  filesystem — irrelevant/must be changed if you're not using the nested-download-tree staging
  path at all (the common case: just use a flat `input.protein_dir`/`input.cds_dir` and delete or
  ignore the `staging:` block).
- **`parse`'s `.parsed.done` marker is vestigial** — it's still produced (a leftover
  from the original Phase 0 scaffold), but nothing downstream reads it anymore; `protein_properties`/
  `cds_properties` read the real parsed tables directly.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Roadmap

- [x] Phase 0 — Input staging (`stage_inputs`, optional)
- [x] Phase 1 — QC + parsing (`qc`, `parse`)
- [x] Phase 2a — Protein properties (composition, pI, GRAVY, aggregation, ...)
- [x] Phase 2b — Intrinsic disorder (metapredict)
- [x] Phase 3 — CDS/codon properties (codonW)
- [x] Phase 4 — Cross-genome summaries, effect sizes, sensitivity analysis
- [x] Phase 5 — Static visualizations (457 figures)
- [x] Phase 6a — Interactive dashboard (Overview, Property Explorer, Species View, Effect Sizes,
      Sensitivity)
- [x] Phase 6b — CDS/codon-level dashboard views
- [x] Phase 6b — PCA/clustering dashboard views
- [x] Phase 6b — Cross-property scatter plots
- [x] Phase 6b — Export/download from the dashboard

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Contact

Farah Cisse — ffcisse@berkeley.edu

Project Link: [https://github.com/ffcisse/genome_pipeline](https://github.com/ffcisse/genome_pipeline)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Acknowledgments

- [JGI Mycocosm / Phycocosm](https://mycocosm.jgi.doe.gov/) — source of the example 9-genome red
  algae dataset and the FASTA download service `download_from_jgi.sh` documents.
- [metapredict](https://github.com/idptools/metapredict) — intrinsic disorder prediction (Phase 2b).
- [codonW](https://anaconda.org/bioconda/codonw) — codon usage / ENC / GC content (Phase 3).
- [Plotly.js](https://plotly.com/javascript/) — the interactive dashboard's plotting library.
- [Snakemake](https://snakemake.readthedocs.io/) — workflow orchestration.
- [Best-README-Template](https://github.com/othneildrew/Best-README-Template) — this README's
  structure.

<p align="right">(<a href="#readme-top">back to top</a>)</p>
