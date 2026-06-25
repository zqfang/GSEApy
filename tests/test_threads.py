"""Regression tests for the per-call rayon thread-pool fix in ``src/lib.rs``.

Background
----------
The Rust entry points used to set the number of worker threads with
``env::set_var("RAYON_NUM_THREADS", threads.to_string())``. That value is read
only once -- when rayon lazily builds its *global* pool on the first parallel
call in the process -- so in a long-lived interpreter (a Python console, a
notebook, a web worker) every call after the first silently ignored its
``threads=`` argument. Mutating the process environment was also a data race
when the Rust backend was driven from multiple Python threads.

The fix builds a dedicated ``rayon::ThreadPool`` per call and runs the work
inside ``pool.install(..)`` (with the GIL released via ``py.detach``).
These tests lock in the user-visible contract that follows from it:

1. ``threads`` only affects *speed*, never results -- identical output for any
   thread count, on repeated calls in the same process.
2. The backend is safe to call concurrently from multiple Python threads and
   still returns results identical to a serial baseline.

The exact thread-count honoring is asserted directly in Rust
(``cargo test --features extension-module thread_pool``); from Python we can
only observe that results stay correct and stable, which is what matters to
callers.
"""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from gseapy.__init__ import gsva, prerank, ssgsea

THREAD_COUNTS = [1, 2, 4]
SEED = 123


def _canon(df):
    """Canonicalize result row order for comparison.

    The Rust layer emits one row per (gene set, sample) in ``HashMap`` iteration
    order, which is randomized per ``HashMap`` instance and therefore varies
    between calls -- independently of the thread count. We only care that the
    *values* are thread-count invariant, so sort by the identifying columns
    before comparing.
    """
    id_cols = [c for c in ("Name", "Term") if c in df.columns]
    return df.sort_values(id_cols).reset_index(drop=True)


def _assert_same(a, b):
    pdt.assert_frame_equal(_canon(a), _canon(b))


# --------------------------------------------------------------------------- #
# Fixtures: small, fully synthetic, network-free inputs.
# --------------------------------------------------------------------------- #
@pytest.fixture
def ranking():
    """A small ranked gene list as a (gene, score) DataFrame."""
    genes = [f"G{i:02d}" for i in range(40)]
    scores = np.linspace(5.0, -5.0, num=len(genes))
    return pd.DataFrame({"gene_name": genes, "rank": scores})


@pytest.fixture
def gene_sets():
    return {
        "SET_A": [f"G{i:02d}" for i in range(0, 12)],
        "SET_B": [f"G{i:02d}" for i in range(14, 26)],
        "SET_C": [f"G{i:02d}" for i in range(28, 40)],
    }


@pytest.fixture
def expression():
    """A genes-by-samples expression matrix (synthetic, deterministic)."""
    rng = np.random.default_rng(0)
    genes = [f"G{i:02d}" for i in range(40)]
    samples = [f"S{j}" for j in range(4)]
    data = rng.normal(loc=0.0, scale=1.0, size=(len(genes), len(samples)))
    return pd.DataFrame(data, index=genes, columns=samples)


def _run_prerank(ranking, gene_sets, threads):
    return prerank(
        rnk=ranking.copy(),
        gene_sets=gene_sets,
        min_size=1,
        max_size=100,
        permutation_num=100,
        threads=threads,
        seed=SEED,
        no_plot=True,
        verbose=False,
        outdir=None,
    ).res2d


# --------------------------------------------------------------------------- #
# prerank
# --------------------------------------------------------------------------- #
def test_prerank_results_independent_of_thread_count(ranking, gene_sets):
    """Same seed, different ``threads`` -> byte-identical results."""
    baseline = _run_prerank(ranking, gene_sets, threads=THREAD_COUNTS[0])
    for threads in THREAD_COUNTS[1:]:
        result = _run_prerank(ranking, gene_sets, threads=threads)
        _assert_same(baseline, result)


def test_prerank_threads_honored_after_first_call(ranking, gene_sets):
    """Core regression for the old global-state bug.

    With the previous ``RAYON_NUM_THREADS`` approach the thread count was frozen
    at the first call for the life of the process. Here we interleave thread
    counts -- including returning to the first value -- and require every call
    to keep producing the canonical result.
    """
    baseline = _run_prerank(ranking, gene_sets, threads=1)
    for threads in [4, 2, 1, 4]:
        result = _run_prerank(ranking, gene_sets, threads=threads)
        _assert_same(baseline, result)


def test_prerank_safe_under_python_thread_concurrency(ranking, gene_sets):
    """Driving the Rust backend from many Python threads must be race-free.

    The old code mutated the process environment on every call, which is a data
    race when several threads call in parallel. The per-call pool removes that
    shared state, and ``py.detach`` lets the calls actually run
    concurrently. We require every concurrent result to match the serial
    baseline and the run to complete without deadlock or crash.
    """
    baseline = _run_prerank(ranking, gene_sets, threads=2)

    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(_run_prerank, ranking, gene_sets, 2) for _ in range(8)]
        results = [f.result() for f in futures]

    for result in results:
        _assert_same(baseline, result)


# --------------------------------------------------------------------------- #
# ssGSEA
# --------------------------------------------------------------------------- #
def test_ssgsea_results_independent_of_thread_count(expression, gene_sets):
    baseline = ssgsea(
        data=expression,
        gene_sets=gene_sets,
        min_size=5,
        max_size=100,
        threads=THREAD_COUNTS[0],
        seed=SEED,
        no_plot=True,
        verbose=False,
        outdir=None,
    ).res2d
    for threads in THREAD_COUNTS[1:]:
        result = ssgsea(
            data=expression,
            gene_sets=gene_sets,
            min_size=5,
            max_size=100,
            threads=threads,
            seed=SEED,
            no_plot=True,
            verbose=False,
            outdir=None,
        ).res2d
        _assert_same(baseline, result)


# --------------------------------------------------------------------------- #
# GSVA
# --------------------------------------------------------------------------- #
def test_gsva_results_independent_of_thread_count(expression, gene_sets):
    baseline = gsva(
        data=expression,
        gene_sets=gene_sets,
        min_size=5,
        max_size=100,
        threads=THREAD_COUNTS[0],
        seed=SEED,
        verbose=False,
        outdir=None,
    ).res2d
    for threads in THREAD_COUNTS[1:]:
        result = gsva(
            data=expression,
            gene_sets=gene_sets,
            min_size=5,
            max_size=100,
            threads=threads,
            seed=SEED,
            verbose=False,
            outdir=None,
        ).res2d
        _assert_same(baseline, result)
