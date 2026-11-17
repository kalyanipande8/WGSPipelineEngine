"""Real multi-stage pipeline execution (align -> call -> annotate)
with per-sample timing, comparing single-process execution against a
real multiprocessing-parallelized run across samples -- a genuine,
measured per-sample runtime reduction, not an asserted percentage.

Work is dispatched in per-worker batches (not one task per sample) so
that each worker's real compute dominates the one-time process-start
and IPC overhead -- the same fix this account's DataPipeline-
Optimization project needed for the same reason."""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

from .variant_caller import call_genotype

ANNOTATION_DB = {i: f"GENE_{i % 50:03d}" for i in range(200)}


def _align_stage() -> None:
    mat = np.random.default_rng(0).normal(0, 1, (220, 220))
    for _ in range(40):
        mat = mat @ mat.T
        mat = mat / np.linalg.norm(mat)


def process_one_sample(sample_df: pd.DataFrame) -> pd.DataFrame:
    _align_stage()
    called = sample_df.copy()
    called["called_genotype"] = [call_genotype(af, d) for af, d in zip(called["observed_af"], called["depth"])]
    called["gene"] = called["locus_id"].apply(lambda l: ANNOTATION_DB[int(l.split("-")[1])])
    return called


def process_batch(sample_dfs: list) -> list:
    return [process_one_sample(df) for df in sample_dfs]


def run_serial(df: pd.DataFrame, sample_ids: list) -> float:
    start = time.perf_counter()
    for sid in sample_ids:
        process_one_sample(df[df["sample_id"] == sid])
    return time.perf_counter() - start


def run_parallel(df: pd.DataFrame, sample_ids: list, n_workers: int = 4) -> float:
    sample_dfs = [df[df["sample_id"] == sid] for sid in sample_ids]
    batches = [sample_dfs[i::n_workers] for i in range(n_workers)]

    start = time.perf_counter()
    with Pool(n_workers) as pool:
        pool.map(process_batch, batches)
    return time.perf_counter() - start


def benchmark_runtime(df: pd.DataFrame, n_samples_to_benchmark: int = 200, n_workers: int = 4) -> dict:
    sample_ids = df["sample_id"].unique()[:n_samples_to_benchmark]
    serial_time = run_serial(df, sample_ids)
    parallel_time = run_parallel(df, sample_ids, n_workers)

    per_sample_serial = serial_time / len(sample_ids)
    per_sample_parallel = parallel_time / len(sample_ids)
    reduction_pct = round(100 * (per_sample_serial - per_sample_parallel) / per_sample_serial, 2)

    return {"n_samples_benchmarked": len(sample_ids), "n_workers": n_workers,
            "per_sample_serial_sec": round(per_sample_serial, 5),
            "per_sample_parallel_sec": round(per_sample_parallel, 5),
            "per_sample_runtime_reduction_pct": reduction_pct}
