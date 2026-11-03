# QSEEG

Code for **Quantile spectral regression and simultaneous inference for repeated EEG recordings**.

**Tianbo Chen and Hanbing Zhu**  
School of Big Data and Statistics, Anhui University, Hefei, China

**Corresponding author:** Hanbing Zhu, [zhuhbecnu@163.com](mailto:zhuhbecnu@163.com)

## Paper overview

Repeated electroencephalographic (EEG) recordings contain dependence both within a segment and across segments from the same person. QSEEG studies covariate associations with rank-based quantile spectral summaries while retaining the person as the independent sampling unit. The rank representation describes temporal ordering relative to each segment's amplitude distribution and complements ordinary spectral analysis of amplitude-sensitive power.

The method combines subject-normalized regression weights, a cluster sandwich covariance estimator, and a studentized subject-level wild bootstrap. It accommodates unequal recording counts, segment-level covariates, and paired within-person contrasts. Confidence bands and global tests jointly account for dependence across the selected frequencies and quantile levels.

The primary target is the regression of finite-window log spectral summaries at a common smoothing resolution. The theory establishes asymptotic normality and bootstrap validity on a fixed frequency-quantile grid as the number of independent subjects increases, allowing segment counts and segment lengths to remain fixed or grow under uniform cluster moment conditions. Additional bias and mean-square error conditions, including an explicit regime for the ratios of subject counts, segment counts, and segment lengths, connect this target to a latent spectral regression.

## Simulation studies and EEG applications

The simulation studies assess simultaneous coverage, precision, and sensitivity to changes in temporal dependence. They show why subject weighting and dependence correction address different aspects of inference. In stochastic-volatility settings, noncentral quantile spectra detect dependence changes that can be missed by median-only and ordinary spectral analyses. Ordinary spectral regression remains competitive in the Gaussian autoregressive settings. Null experiments with 40 subjects give rejection probabilities of 0.064 in the regression-surface model and 0.050 in the stochastic-volatility model. Surface-model coverage is 0.936, slightly below the nominal 0.95.

Applications use the Sleep-EDF Database Expanded. The age study examines adjusted between-person associations, including quantile-dependent age associations in parietal-occipital recordings. The paired study examines frequency-dependent temazepam-minus-placebo contrasts within the same person and sleep stage. These analyses describe spectral associations and paired contrasts, rather than clinical treatment benefits.

## Download and reproduce

The computational materials are distributed in two parts:

1. Download **QSEEG_code.zip** from the [data bundle release](https://github.com/tianbochen1/QSEEG/releases/tag/v1.0.1). It contains source metadata and sleep-stage annotations, derived spectral features for 197 recordings, and the existing simulation and EEG results.
2. Download the [current repository source archive](https://github.com/tianbochen1/QSEEG/archive/refs/heads/main.zip). It contains the current programs and manuscript components, together with additional null-simulation records.

Extract both archives. Merge the repository's `research/` directory into `QSEEG_code/research/`, replacing matching files and retaining the other data and result files from the data bundle. Then enter the `QSEEG_code` directory and run:

```bash
python -m pip install -r requirements.txt
python research/code/report_simulations.py
python research/code/report_eeg.py
python research/code/build_manuscript.py
```

These commands reproduce the tables and figures and assemble the current LaTeX manuscript. Python 3.12 is used for the analyses; pdfLaTeX is needed to compile the manuscript. See [reproduction instructions](docs/REPRODUCIBILITY.md) for recalculating the simulations, EEG regressions, and spectral features.

The GitHub-generated **Source code** archives contain the repository files only. Both downloads are needed to reproduce the current manuscript from the supplied numerical results.

## Repository contents

| Location | Contents |
| --- | --- |
| `research/code/` | Twelve Python programs for spectral features, estimation, inference, simulations, EEG analysis, and reporting |
| `research/data/` | Official participant metadata, recording list, and source-file checksum manifest |
| `research/manuscript/` | Manuscript components and bibliography |
| `research/results/` | Additional null-simulation records |
| `requirements.txt` | Python dependencies |
| `docs/REPRODUCIBILITY.md` | Data provenance and reproduction instructions |

The larger derived data and result files are distributed in the release bundle.

## Data source

The analyses use [Sleep-EDF Database Expanded, version 1.0.0](https://physionet.org/content/sleep-edfx/1.0.0/), DOI [10.13026/C2X676](https://doi.org/10.13026/C2X676). The source data are distributed under the Open Data Commons Attribution License v1.0. Data citations and acquisition details are given in the reproduction instructions.


