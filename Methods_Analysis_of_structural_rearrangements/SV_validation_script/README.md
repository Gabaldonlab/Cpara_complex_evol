# General information

Coverage plots and a coverage-shift test for structural variants (SVs) that differ
between a target sample and background samples.

For each SV the script:

1. Builds windows around the breakpoints (one block for intra-chromosomal SVs,
   one block per breakend for inter-chromosomal SVs), each split into a fixed number of bins.
2. Computes mean depth per bin with mosdepth and divides by each sample's
   genome-wide mean depth, so 1 = genome-wide average coverage.
3. Plots all samples (target first) with breakpoints marked.
4. For real intra-chromosomal SVs only, tests whether the target-vs-background
   log2 coverage ratio inside the event differs from the flanks.

## Requirements

Tested on Linux with:

| Software | Version | Notes |
|---|---|---|
| Python | 3.10.19 | ≥ 3.10 required |
| mosdepth | 0.3.11 | must be on `PATH` |
| numpy | 2.2.6 | |
| pandas | 2.3.3 | |
| matplotlib | 3.10.7 | |
| reportlab | 4.4.4 | optional, for the PDF generation|
| svglib | 1.6.0 | optional, for the PDF generation |

To create a matching environment with conda:

```bash
conda env create -f SVanalysis_env.yml
conda activate SVanalysis_env
```

## Usage

```bash
python SVs_validation.py \
    --sv_table SVs_changing.tab \
    --paths_table paths.tsv \
    --outdir sv_coverage_out \
    --threads 4
```

Run `python SVs_validation.py --help` for all options and defaults.

## Inputs

**`--sv_table`** (tab-separated), e.g. `SVs_changing.tab` from `integrate_several_samples`:

| Column | Required | Description |
|---|---|---|
| `INFO_variantID` | yes | SV ID, `TYPE\|chrA:posA-chrB:posB` or `TYPE\|chr:start-end` |
| `target_sampleID` | yes | Sample in which the SV is differential |
| `type_background` | yes | Only rows with `low_confidence` are used |
| `background_sampleIDs` | no | Comma-separated background samples |
| `type_differetial_variant` | no | Label used in file names (spelling matches the upstream column) |

Contig names must not contain `:` or `|`. SV IDs containing `like`
(case-insensitive) are treated as "like" SVs: plotted, but not included in the statistical test. Only "real" SVs are tested.

By default only the first row per `INFO_variantID` is kept. Use `--all_targets`
to process every (SV, target) pair.

**`--paths_table`** (tab-separated): columns `sampleID`, `sorted_bam`
(coordinate-sorted, indexed BAM). Sample IDs must be unique.

## Method

**Windows.** Intra-chromosomal SVs get one block from `start − flank` to
`end + flank`, where flank = max(`intra_flank_pct` × SV length, `intra_min_flank_bp`),
optionally capped by `intra_max_flank_bp`. Inter-chromosomal SVs get one block of
± `inter_flank_bp` around each breakend. Blocks are clamped to contig ends and split
into `bins_per_block` equal bins.

**Normalization.** Bin depth is divided by the sample's genome-wide mean depth.
The meaning of a given value depends on ploidy: in a haploid genome,
2 indicates one extra copy; in a diploid, 1.5.

**Test.** For each bin, log2(target / median of backgrounds) is computed; with no
background samples, the target is compared with 1. The effect is
median(event bins) − median(flank bins). Flank bins are taken next to the event,
skipping `flank_gap_bins` on each side, up to as many as there are event bins.
Significance is calculated with exact two-sided test over all `n − 1` non-identity
circular shifts of the bin series, with event and flank masks held fixed:

    p = (number of shifts with |effect| ≥ |observed| + 1) / n

where n = number of bins in the block. The smallest achievable p is 1/n
(0.02 with the default 50 bins) and is reported as `p_min_achievable`.
Benjamini–Hochberg adjustment is applied across all tested blocks in
`stats/coverage.stats.tsv` only.

## Outputs

| Path | Content |
|---|---|
| `beds/<sv>.windows.bed` | Windows (name = `b<block>w<bin>`) |
| `mosdepth/<sample>.mosdepth.summary.txt` | Genome-wide summary (cached) |
| `mosdepth/<sv>__<sample>.regions.bed.gz` | Per-window depth (cached) |
| `mosdepth/*.cache.json` | Cache stamps (BAM path, size, mtime) |
| `figs/real/`, `figs/like/` | `<sv>_<target>_<diff>.coverage.svg` |
| `figs/coverage_plots_REAL_SVs.pdf` | Real SVs, 2 per A4 page: tested ones by adjusted p, then untested |
| `stats/coverage_stats/<sv>_<target>_<diff>.coverage.stats.tsv` | Per-SV results |
| `stats/coverage_stats/by_sample/<sample>.coverage.stats.tsv` | Per-target results sorted by raw p (no BH) |
| `stats/permutations/<sv>_<target>_<diff>.permutations.tsv` | Observed effect (`shift` = 0) and all shifted effects |
| `stats/coverage.stats.tsv` | All tested blocks with BH-adjusted p |

Columns of `stats/coverage.stats.tsv`:

| Column | Description |
|---|---|
| `effect_median_diff_log2` | median(event) − median(flanks) of log2(target/background) |
| `p_perm_two_sided` | Exact circular-shift p-value |
| `p_min_achievable` | 1 / number of bins |
| `p_adj_bh` | Benjamini–Hochberg adjusted p across all tested blocks |
| `n_windows`, `n_event_bins`, `n_flank_bins` | Bin counts used |
| `figure`, `cov_stats_file` | Paths relative to `--outdir` |