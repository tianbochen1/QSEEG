"""Null calibration using the covariates and residual covariance of the paired EEG design."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import argparse
import json
import numpy as np
import pandas as pd
from analyse_eeg import load
from core import fit

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'paired_calibration'
SETTINGS = [('Gaussian', False), ('Gaussian', True), ('t5', False), ('t5', True)]


def design():
    df, q, _, _, _ = load('ST', 2, 0, 'main')
    complete = df.groupby('subject').treatment.nunique()
    rows, responses, ages = [], [], []
    for subject in complete[complete == 2].index:
        exposed = (df.subject == subject) & (df.treatment == 1)
        control = (df.subject == subject) & (df.treatment == 0)
        rows.append([1., df.loc[exposed, 'night'].iloc[0] - df.loc[control, 'night'].iloc[0],
                     df.loc[exposed, 'elapsed'].mean() - df.loc[control, 'elapsed'].mean()])
        responses.append((q[exposed].mean(0) - q[control].mean(0)).ravel())
        ages.append(df.loc[exposed, 'age'].iloc[0])
    X, Y, age = np.asarray(rows), np.asarray(responses), np.asarray(ages)
    beta = np.linalg.solve(X.T @ X, X.T @ Y)
    residual = Y - X @ beta
    factor = residual / np.sqrt(len(X) - X.shape[1])
    mean = X[:, 1:] @ beta[1:]
    scale = np.exp(.35 * (age - age.mean()) / age.std(ddof=0))
    scale /= np.sqrt(np.mean(scale ** 2))
    return X, mean, factor, scale


def run_setting(args):
    setting, repetitions, bootstrap = args
    distribution, heterogeneous = SETTINGS[setting]
    X, mean, factor, scale = design()
    N = len(X)
    rng = np.random.default_rng(731426 + 100003 * setting)
    path = OUT / f'setting_{setting}.jsonl'
    results = []
    with path.open('w') as stream:
        for replication in range(repetitions):
            error = rng.standard_normal((N, N)) @ factor
            if distribution == 't5':
                error *= np.sqrt(3 / rng.chisquare(5, size=N))[:, None]
            if heterogeneous:
                error *= scale[:, None]
            fit_result = fit(X, mean + error, np.arange(N), coef=0,
                             bootstrap=bootstrap, rng=rng)
            reject = fit_result['p'] < .05
            cover = np.all(np.abs(fit_result['effect']) <= fit_result['critical'] * fit_result['se'])
            row = dict(replication=replication, p=fit_result['p'], reject=bool(reject),
                       cover=bool(cover), width=float(np.mean(2 * fit_result['critical'] * fit_result['se'])))
            stream.write(json.dumps(row) + '\n')
            results.append(row)
    return dict(distribution=distribution, heterogeneous=heterogeneous, N=N,
                repetitions=repetitions, bootstrap=bootstrap,
                rejection=np.mean([r['reject'] for r in results]),
                coverage=np.mean([r['cover'] for r in results]),
                width=np.mean([r['width'] for r in results]))


def report():
    summaries = []
    for setting, (distribution, heterogeneous) in enumerate(SETTINGS):
        rows = [json.loads(line) for line in (OUT / f'setting_{setting}.jsonl').read_text().splitlines()]
        summaries.append(dict(distribution=distribution, heterogeneous=heterogeneous,
                              repetitions=len(rows), rejection=np.mean([r['reject'] for r in rows]),
                              coverage=np.mean([r['cover'] for r in rows]),
                              width=np.mean([r['width'] for r in rows])))
    pd.DataFrame(summaries).to_csv(OUT / 'summary.csv', index=False)
    lines = [r'\begin{table}[htbp]', r'\centering',
             r'\caption{Null calibration for the 22-subject paired design. The covariates and grid covariance are based on the stage-2 Fpz--Cz analysis. Errors have a Gaussian or standardized multivariate $t_5$ distribution, with a common or age-dependent subject scale. Rejection and simultaneous coverage are evaluated over the full frequency--quantile grid.}\label{tab:paired-calibration}',
             r'\begin{tabular}{llrr}', r'\toprule',
             r'Error distribution & Subject scale & Rejection & Coverage \\', r'\midrule']
    for row in summaries:
        label = 'Gaussian' if row['distribution'] == 'Gaussian' else '$t_5$'
        scale = 'Age-dependent' if row['heterogeneous'] else 'Common'
        lines.append(f"{label} & {scale} & {row['rejection']:.3f} & {row['coverage']:.3f} " + r'\\')
    lines += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
    (ROOT / 'manuscript' / 'table_paired_calibration.tex').write_text('\n'.join(lines) + '\n')
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--repetitions', type=int, default=1000)
    parser.add_argument('--bootstrap', type=int, default=399)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--report-only', action='store_true')
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if not args.report_only:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for result in pool.map(run_setting, [(s, args.repetitions, args.bootstrap) for s in range(4)]):
                print(json.dumps(result), flush=True)
        (OUT / 'settings.json').write_text(json.dumps(dict(N=22, repetitions=args.repetitions,
                                                          bootstrap=args.bootstrap, age_scale_coefficient=.35), indent=2) + '\n')
    report()
