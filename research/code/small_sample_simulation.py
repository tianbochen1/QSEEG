"""Null calibration with 20 and 40 independent subjects."""

import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import argparse
import concurrent.futures as cf
import csv
import json
import math
from pathlib import Path
import time

from simulate import one


ROOT = Path(__file__).resolve().parents[1]
CONFIGS = [
    dict(kind="surface", N=20, m=4, effect=0, id=108),
    dict(kind="sv", N=20, m=6, n=512, effect=0, id=109),
    dict(kind="surface", N=40, m=4, effect=0, id=110),
    dict(kind="sv", N=40, m=6, n=512, effect=0, id=111),
]
REPLICATIONS = 1000
METHODS = {
    "surface": ("Cluster", "Average", "Independent"),
    "sv": ("Quantile", "Median", "Ordinary", "Independent"),
}


def read_results(path):
    records = {}
    configs = {config["id"]: config for config in CONFIGS}
    if path.exists():
        for line in path.read_text().splitlines():
            record = json.loads(line)
            config = record["config"]
            assert config == configs[config["id"]], config
            key = (config["id"], record["rep"])
            assert 0 <= record["rep"] < REPLICATIONS, key
            assert key not in records, key
            assert set(record["metrics"]) == set(METHODS[config["kind"]])
            for metric in record["metrics"].values():
                assert metric["reject"] in (0, 1)
                assert all(math.isfinite(value) for value in metric.values())
                if config["kind"] == "surface":
                    assert metric["cover"] == 1 - metric["reject"]
            records[key] = record
    return records


def wilson_interval(successes, total):
    z = 1.959963984540054
    probability = successes / total
    denominator = 1 + z * z / total
    center = (probability + z * z / (2 * total)) / denominator
    half = z * math.sqrt(
        probability * (1 - probability) / total + z * z / (4 * total * total)
    ) / denominator
    return max(0, center - half), min(1, center + half)


def summarize(records, destination):
    expected = {
        (config["id"], rep)
        for config in CONFIGS
        for rep in range(REPLICATIONS)
    }
    assert set(records) == expected
    rows = []
    for config in CONFIGS:
        for method in METHODS[config["kind"]]:
            metrics = [
                records[(config["id"], rep)]["metrics"][method]
                for rep in range(REPLICATIONS)
            ]
            rejects = sum(metric["reject"] for metric in metrics)
            lower, upper = wilson_interval(rejects, REPLICATIONS)
            row = dict(
                **config,
                method=method,
                replications=REPLICATIONS,
                bootstrap_draws=399,
                rejection_count=rejects,
                reject=rejects / REPLICATIONS,
                reject_mcse=math.sqrt(
                    (rejects / REPLICATIONS) * (1 - rejects / REPLICATIONS)
                    / REPLICATIONS
                ),
                reject_wilson_lower=lower,
                reject_wilson_upper=upper,
                width=sum(metric["width"] for metric in metrics) / REPLICATIONS,
            )
            if config["kind"] == "surface":
                covers = sum(metric["cover"] for metric in metrics)
                lower, upper = wilson_interval(covers, REPLICATIONS)
                row.update(
                    coverage_count=covers,
                    cover=covers / REPLICATIONS,
                    cover_wilson_lower=lower,
                    cover_wilson_upper=upper,
                    rmse=math.sqrt(
                        sum(metric["rmse"] ** 2 for metric in metrics)
                        / REPLICATIONS
                    ),
                )
            rows.append(row)
    fields = list(dict.fromkeys(field for row in rows for field in row))
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    destination = ROOT / "results" / "small_sample_final.jsonl"
    destination.parent.mkdir(exist_ok=True)
    records = read_results(destination)
    jobs = [
        (config, rep)
        for config in CONFIGS
        for rep in range(REPLICATIONS)
        if (config["id"], rep) not in records
    ]
    start = time.monotonic()
    print(f"Pending replications: {len(jobs)}", flush=True)
    with cf.ProcessPoolExecutor(max_workers=args.workers) as executor:
        with destination.open("a") as output:
            for completed, record in enumerate(executor.map(one, jobs, chunksize=2), 1):
                output.write(json.dumps(record) + "\n")
                output.flush()
                if completed % 100 == 0:
                    print(
                        f"Completed {completed}; elapsed {time.monotonic() - start:.1f} s",
                        flush=True,
                    )
    rows = summarize(
        read_results(destination), ROOT / "results" / "small_sample_summary.csv"
    )
    print(json.dumps(rows, indent=2), flush=True)
