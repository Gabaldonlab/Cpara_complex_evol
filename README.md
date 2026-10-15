**General information**

This repository contains the code and resources for the analysis of structural rearrangements, copy number variants and aneuploidies described in our paper

> **Reversible haploidisation and convergent genomic routes to antifungal resistance in the *Candida parapsilosis* species complex**
>
> Juan Carlos Nunez-Rodriguez†, Vladislav Biriukov†, Álvaro Redondo-Río†, Mrinalini Parmar, Timothy James Nickels, Miquel Àngel Schikora-Tamarit, Islam Ahaik, Valentina del Olmo, Toni Gabaldón\*
>
> † These authors contributed equally. \* Corresponding author.

## Repository structure

```
Cpara_complex_evol/
├── README.md
├── Methods_Analysis_of_structural_rearrangements/
│   ├── SV_parameters/                        # SV calling parameter files per species
│   │   ├── <species>_sv_calling_parameters.json          # selected configuration (used in the study)
│   │   ├── <species>_sv_calling_parameters_strict.json   # stringent configuration (benchmarked, rejected)
│   │   └── SV_analyze_parameters_heatmaps/   # benchmark results on simulated genomes
│   ├── repeats/                              # filtered repeat annotations per species
│   └── SV_validation_script/                 # coverage-based curation of differential SVs
└── Methods_Detection_of_copy_number_variants_and_aneuploidies/
    └── <sample>_aneuploidies.bed             # manually curated aneuploidies per sample
```

## Contents

### `SV_parameters/`

Three configurations were benchmarked against truth sets of
rearrangements simulated on each reference genome:

1. **Defaults:** GRIDSS/perSVade default parameters.
2. **Selected:** `<species>_sv_calling_parameters.json`, used for all analyses in the paper.
3. **Stringent:** `<species>_sv_calling_parameters_strict.json`, deposited for reproducibility, but rejected on the benchmark.

The selected and stringent configurations are identical in all 21 parameters except
`wrong_FILTERtags`, the set of GRIDSS FILTER tags used to discard calls:

| Configuration | Calls discarded if tagged with |
|---|---|
| Selected | `LOW_QUAL`, `INSUFFICIENT_SUPPORT` |
| Stringent | the above, plus `ASSEMBLY_ONLY`, `NO_ASSEMBLY`, `ASSEMBLY_TOO_FEW_READ`, `ASSEMBLY_TOO_SHORT`, `REF`, `SINGLE_ASSEMBLY` |

#### `SV_parameters/SV_analyze_parameters_heatmaps/`

Heatmaps of the benchmark results for each species of the SV parameters. 

### `repeats/`

Final repeat annotations used for SV calling, one per reference genome
(`<species>_simple_repeats.tab`). Repeats were inferred with RepeatMasker through the
`infer_repeats` module of [perSVade](https://github.com/Gabaldonlab/perSVade). The file was then filtered  and only entries whose RepeatMasker class is
`Simple_repeat` or `Low_complexity` were kept.

### `SV_validation_script/`

`SVs_validation.py` is the custom curation script referred to in the paper. It
generates read-depth coverage plots around differential SVs obtained by
[perSVade](https://github.com/Gabaldonlab/perSVade) and tests intra-chromosomal SVs
for a coverage shift between the event and its flanks, using an exact circular-shift
permutation test. These plots and statistics were used to identify possible artefacts
among differential SV calls. See
[`SV_validation_script/README.md`](Methods_Analysis_of_structural_rearrangements/SV_validation_script/README.md)
for requirements, usage, method and outputs. The conda environment is provided in
`SVanalysis_env.yaml`.

### `Methods_Detection_of_copy_number_variants_and_aneuploidies/`

Aneuploidies defined for each sample (`<sample>_aneuploidies.bed`), one file per
sample in which aneuploidies were detected.

## Data availability

Raw reads from this study are deposited in the NCBI Sequence Read Archive under BioProject **PRJNA1280132**.