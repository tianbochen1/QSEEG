# QSEEG

Code for **Quantile spectral regression and simultaneous inference for repeated EEG recordings**.

**Tianbo Chen and Hanbing Zhu**  
School of Big Data and Statistics, Anhui University, Hefei, China

**Corresponding author:** Hanbing Zhu, [zhuhbecnu@163.com](mailto:zhuhbecnu@163.com)

## Paper overview

Repeated electroencephalographic (EEG) recordings contain dependence both within a segment and across segments from the same person. QSEEG studies covariate associations with rank-based quantile spectral summaries while retaining the person as the independent sampling unit. The rank representation describes temporal ordering relative to each segment's amplitude distribution and complements ordinary spectral analysis of amplitude-sensitive power.

The method combines subject-normalized regression weights, a cluster sandwich covariance estimator, and a studentized subject-level wild bootstrap. At each channel and window length, the analyses retain segments with positive rank-indicator variance at every tested quantile. The regression accommodates unequal recording counts, segment-level covariates, and paired within-person contrasts. Confidence bands and global tests jointly account for dependence across the selected frequencies and quantile levels.

The finite-resolution regression target is defined by the segment length, spectral smoothing, and analysis grid. The theory establishes asymptotic normality and bootstrap validity on a fixed frequency-quantile grid as the number of independent subjects increases, allowing segment counts and segment lengths to remain fixed or grow under uniform cluster moment conditions. Additional conditional bias and mean-square error assumptions give sufficient growth regimes linking this target to a latent spectral regression; these rates require verification for the chosen time-series model.

## Simulation studies and EEG applications

The simulation studies assess simultaneous coverage, precision, and sensitivity to changes in temporal dependence. They show why subject weighting and dependence correction address different aspects of inference. In stochastic-volatility settings, nonmedian quantile spectra detect dependence changes that can be missed by median-only and ordinary spectral analyses. Ordinary spectral regression remains competitive in the Gaussian autoregressive settings. Null experiments with 40 subjects give rejection probabilities of 0.064 in the regression-surface model and 0.050 in the stochastic-volatility model. The corresponding simultaneous coverages are 0.936 and 0.950. Calibration experiments based on the 22-subject paired EEG design give rejection probabilities of 0.038–0.053 and simultaneous coverage of 0.947–0.962 across four Gaussian and heavy-tailed error settings.

Applications use the Sleep-EDF Database Expanded. Within the parietal-occipital lead, the estimated age association varies across quantile levels. The paired study examines frequency-dependent temazepam-minus-placebo contrasts within the same person and sleep stage. The paired study detected an overall contrast, with insufficient evidence of variation across quantile levels after adjustment.

## Download and reproduce

Download **[QSEEG_code.zip from version 1.0.5](https://github.com/tianbochen1/QSEEG/releases/download/v1.0.5/QSEEG_code.zip)**. This complete bundle contains the matching programs and manuscript components, source metadata and annotations, derived spectral features for 197 recordings, and all numerical results needed to reproduce the paper. No merging of separate source and data archives is required.

Extract the bundle, enter `QSEEG_code`, and run:

```bash
python -m pip install -r requirements.txt
python research/code/report_simulations.py
python research/code/report_eeg.py
python research/code/paired_calibration.py --report-only
python research/code/build_manuscript.py
```

These commands reproduce the tables and figures and assemble the LaTeX manuscript. Python 3.12 is used for the analyses; pdfLaTeX is needed to compile the manuscript. See [reproduction instructions](docs/REPRODUCIBILITY.md) for recalculating the simulations, EEG regressions, and spectral features.

Use the **QSEEG_code.zip** release asset for the complete materials. GitHub's automatically generated source archives contain only repository files. The release includes a SHA256 checksum, and [the content manifest](docs/code_manifest.json) records the checksum of every file in the computational bundle.

## Repository contents

| Location | Contents |
| --- | --- |
| `research/code/` | Thirteen Python programs for spectral features, estimation, inference, simulations, EEG analysis, and reporting |
| `research/data/` | Official participant metadata, recording list, and source-file checksum manifest |
| `research/manuscript/` | Manuscript components and bibliography |
| `research/results/` | Additional null-simulation and paired-calibration records |
| `requirements.txt` | Python dependencies |
| `docs/REPRODUCIBILITY.md` | Data provenance and reproduction instructions |

The larger derived data and result files are distributed in the release bundle.

## Data source

The analyses use [Sleep-EDF Database Expanded, version 1.0.0](https://physionet.org/content/sleep-edfx/1.0.0/), DOI [10.13026/C2X676](https://doi.org/10.13026/C2X676). The source data are distributed under the Open Data Commons Attribution License v1.0. Data citations and acquisition details are given in the reproduction instructions.



