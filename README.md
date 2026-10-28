# QSEEG

Code for **Quantile spectral regression and simultaneous inference for repeated EEG recordings**.

**Tianbo Chen and Hanbing Zhu**  
School of Big Data and Statistics, Anhui University, Hefei, China

**Corresponding author:** Hanbing Zhu, [zhuhbecnu@163.com](mailto:zhuhbecnu@163.com)

## Paper overview

Repeated electroencephalographic (EEG) recordings contain dependence both within a segment and across segments from the same person. QSEEG studies covariate associations with rank-based quantile spectral summaries while retaining the person as the independent sampling unit. The rank representation describes temporal ordering relative to each segment's amplitude distribution and complements ordinary spectral analysis of amplitude-sensitive power.

The method combines subject-normalized regression weights, a cluster sandwich covariance estimator, and a studentized subject-level wild bootstrap. It accommodates unequal recording counts, segment-level covariates, and paired within-person contrasts. Confidence bands and global tests jointly account for dependence across the selected frequencies and quantile levels.

The primary target is the regression of finite-window log spectral summaries at a common smoothing resolution. The theory establishes asymptotic normality and bootstrap validity on a fixed frequency-quantile grid as the number of independent subjects increases. Additional approximation conditions connect this target to a latent spectral regression.

## Simulation studies and EEG applications

The simulation studies assess simultaneous coverage, precision, and sensitivity to changes in temporal dependence. They show why subject weighting and dependence correction address different aspects of inference. In stochastic-volatility settings, noncentral quantile spectra detect dependence changes that can be missed by median-only and ordinary spectral analyses. Ordinary spectral regression remains competitive in the Gaussian autoregressive settings.

Applications use the Sleep-EDF Database Expanded. The age study examines adjusted between-person associations, including quantile-dependent age associations in parietal-occipital recordings. The paired study examines frequency-dependent temazepam-minus-placebo contrasts within the same person and sleep stage. These analyses describe spectral associations and paired contrasts, rather than clinical treatment benefits.

## Download and reproduce

Download **QSEEG_code.zip** from the [latest release](https://github.com/tianbochen1/QSEEG/releases/latest). The complete bundle contains the analysis programs, source metadata and sleep-stage annotations, derived spectral features for 197 recordings, numerical results, manuscript source components, and English reproduction instructions.

After extraction, enter the `QSEEG_code` directory and run:

```bash
python -m pip install -r requirements.txt
python research/code/report_simulations.py
python research/code/report_eeg.py
python research/code/build_manuscript.py
```

These commands reproduce the tables and figures from the included numerical results and assemble the LaTeX manuscript. Python 3.12 is used for the analyses; pdfLaTeX is needed to compile the manuscript. See [reproduction instructions](docs/REPRODUCIBILITY.md) for recalculating the simulations, EEG regressions, and spectral features.

The GitHub-generated **Source code** archives contain the repository files only. Use **QSEEG_code.zip** for the complete bundle with derived data and results.

## Repository contents

| Location | Contents |
| --- | --- |
| `research/code/` | Eleven Python programs for spectral features, estimation, inference, simulations, EEG analysis, and reporting |
| `research/data/` | Official participant metadata, recording list, and source-file checksum manifest |
| `research/manuscript/` | Manuscript components and bibliography |
| `requirements.txt` | Python dependencies |
| `docs/REPRODUCIBILITY.md` | Data provenance and reproduction instructions |

The larger derived data and result files are distributed in the release bundle.

## Data source

The analyses use [Sleep-EDF Database Expanded, version 1.0.0](https://physionet.org/content/sleep-edfx/1.0.0/), DOI [10.13026/C2X676](https://doi.org/10.13026/C2X676). The source data are distributed under the Open Data Commons Attribution License v1.0. Data citations and acquisition details are given in the reproduction instructions.
