# QSEEG

Computational materials for *Quantile spectral regression and simultaneous inference for repeated EEG recordings*.

**Authors:** Tianbo Chen and Hanbing Zhu. Both authors are affiliated with the School of Big Data and Statistics, Anhui University, Hefei, China.

**Corresponding author:** Hanbing Zhu, zhuhbecnu@163.com.

**Repository:** https://github.com/tianbochen1/QSEEG

## Software

The analyses use Python 3.12 with NumPy, SciPy, pandas, Matplotlib, pyedflib, and xlrd. Install the dependencies from the top-level directory of this archive:

```bash
python -m pip install -r requirements.txt
```

A standard TeX Live installation with pdfLaTeX is needed only to compile the manuscript. The bibliography is embedded in the assembled LaTeX source in apalike format.

## Reproduce the tables, figures, and manuscript

The archive includes the numerical simulation results, EEG regression results, and spectral features for all 197 recordings. Run the following commands from the directory containing `research/`:

```bash
python research/code/report_simulations.py
python research/code/report_eeg.py
python research/code/build_manuscript.py
```

The assembled source is `output/QSEEG/QSEEG.tex`, with seven figures in `output/QSEEG/frame/`. Compile from that directory:

```bash
cd output/QSEEG
pdflatex -interaction=nonstopmode -halt-on-error QSEEG.tex
pdflatex -interaction=nonstopmode -halt-on-error QSEEG.tex
pdflatex -interaction=nonstopmode -halt-on-error QSEEG.tex
```

The separately supplied manuscript archive contains only this single LaTeX file and its figure directory.

## Recalculate the analyses

The simulation and EEG scripts reuse completed results. To recalculate all analyses, first preserve the supplied results under a different directory name. From the top-level directory, run:

```bash
mv research/results research/reference_results
mkdir research/results
python research/code/simulate.py --reps 1000 --workers 5
python research/code/extra_simulation.py
python research/code/gaussian_comparison.py
python research/code/analyse_eeg.py
python research/code/report_simulations.py
python research/code/report_eeg.py
python research/code/build_manuscript.py
```

The supplied feature archives are sufficient for recalculating the EEG regressions. The Gaussian bootstrap comparison replaces its own output when run. The other simulation scripts resume incomplete calculations using configuration and replication identifiers. Random-number generators are initialized in the analysis scripts so the calculations can be reproduced. Minor floating-point differences can occur across numerical-library versions and platforms.

To recalculate only the EEG analyses, move `research/results/eeg/` to a separate directory, then run `analyse_eeg.py` and `report_eeg.py` using the commands above.

## Recreate spectral features from source recordings

Source data are from the Sleep-EDF Database Expanded, version 1.0.0:

- Data: https://physionet.org/content/sleep-edfx/1.0.0/
- DOI: https://doi.org/10.13026/C2X676
- Data license: Open Data Commons Attribution License v1.0, https://opendatacommons.org/licenses/by/1-0/

The archive includes the official recording list, SHA256 manifest, participant spreadsheets, sleep-stage annotations, and derived spectral features. Full polysomnographic signals are obtained from the public source when needed. To rebuild the spectral features, first move the supplied `research/data/features/` directory to a separate location, then run:

```bash
python research/code/acquire_eeg.py
python research/code/eeg_features.py
```

Acquisition verifies the raw recording and annotation hashes against the source manifest. It selects nonoverlapping, stage-specific 60-second intervals by their positions in the annotated recording, with up to 20 intervals per stage and night. Both Fpz-Cz and Pz-Oz channels are retained. Channel-specific quality indicators are stored with the selected intervals. Each raw recording is removed from the acquisition directory after its selected intervals have been extracted; the extracted intervals remain available for feature calculation.

The annotation reader uses the standard EDF date and time fields and the annotation stream. This accommodates seven source annotations with inconsistent dates in the additional recording-description field while preserving the original annotation files.

Please cite the database and the associated papers:

- Mourtazaev, M. S., Kemp, B., Zwinderman, A. H., and Kamphuisen, H. A. C. (1995). Age and gender affect different characteristics of slow waves in the sleep EEG. *Sleep*, 18, 557-564. https://doi.org/10.1093/sleep/18.7.557
- Kemp, B., Zwinderman, A. H., Tuk, B., Kamphuisen, H. A. C., and Oberye, J. J. L. (2000). Analysis of a sleep-dependent neuronal feedback loop: The slow-wave microcontinuity of the EEG. *IEEE Transactions on Biomedical Engineering*, 47, 1185-1194. https://doi.org/10.1109/10.867928
- Goldberger, A. L., et al. (2000). PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals. *Circulation*, 101, e215-e220. https://doi.org/10.1161/01.CIR.101.23.e215

## Analysis conventions

The primary response is a finite-window, smoothed log rank spectrum on a common frequency-quantile grid. Subject-normalized weights determine the regression target. Cluster sandwich standard errors and studentized subject-level wild bootstrap draws retain dependence across segments, nights, frequencies, and quantile levels.

The age analysis uses adjusted between-person associations. The medication analysis uses within-person temazepam-minus-placebo contrasts. Stage-2 Fpz-Cz is primary in both studies. REM and Pz-Oz analyses are secondary. Direct quantile heterogeneity tests compare all nonmedian quantile curves with the median curve, with Holm adjustment across the four stage-channel combinations within each study.

There are 25 main and additional simulation configurations with 1,000 datasets per configuration. The manuscript presents 24 configurations. The archived contamination alternative is additional material without a matched null experiment. Four further configurations compare a Gaussian multiplier construction with fixed standard errors against the studentized wild bootstrap. All simulation summaries use complete replication sets. Root mean squared error is calculated by averaging squared grid errors across replications before taking the square root. Monte Carlo probability error bars use Wilson intervals.

## Files

| Location | Contents |
| --- | --- |
| `research/code/core.py` | Spectral responses, weighted regression, covariance estimation, and bootstrap inference |
| `research/code/simulate.py` | Main regression-surface and raw time-series simulation designs |
| `research/code/extra_simulation.py` | Additional effect sizes and small-subject designs |
| `research/code/gaussian_comparison.py` | Gaussian multiplier comparison |
| `research/code/read_annotations.py` | Sleep-stage annotation reader |
| `research/code/acquire_eeg.py` | Recording acquisition, interval selection, and quality indicators |
| `research/code/eeg_features.py` | Spectral feature extraction at the main and sensitivity resolutions |
| `research/code/analyse_eeg.py` | Age, paired medication, and resolution sensitivity analyses |
| `research/code/report_simulations.py` | Simulation tables and figures |
| `research/code/report_eeg.py` | EEG tables and figures |
| `research/code/build_manuscript.py` | Assembly of the single-file manuscript |
| `research/data/` | Source metadata, annotations, and 197 derived feature archives |
| `research/results/` | Numerical results and inclusion tables |
| `research/manuscript/` | English manuscript components and bibliography |

Bootstrap critical values use the same strict rejection convention as the Monte Carlo tail probabilities. The bootstrap variance calculation is algebraically equivalent to refitting each resampled weighted regression.
