import os
import re
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXE = ROOT / "build" / "mpich-release" / "OpenSeesMP"
BENCHMARK = ROOT / "EXAMPLES" / "ParallelModelMP" / "metisPartitionBenchmark.tcl"
ROW = re.compile(
    r"METIS_BENCHMARK rank=(?P<rank>\d+) ranks=(?P<ranks>\d+) "
    r"before=(?P<before>\d+)\s+after=(?P<after>\d+)"
)


def executable():
    path = Path(os.environ.get("OPENSEESMP", DEFAULT_EXE))
    if not path.is_file():
        pytest.skip(f"OpenSeesMP executable not found: {path}")
    return path


def run_benchmark(ranks, objective="volume", seed=42):
    env = os.environ.copy()
    env.setdefault("TCL_LIBRARY", "/opt/homebrew/opt/tcl-tk@8/lib/tcl8.6")
    result = subprocess.run(
        ["mpiexec", "-n", str(ranks), str(executable()), str(BENCHMARK),
         objective, str(seed)],
        cwd=ROOT,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=30,
    )
    return result, [
        {key: int(value) for key, value in match.groupdict().items()}
        for match in ROW.finditer(result.stdout)
    ]


@pytest.mark.parametrize("ranks", [2, 4, 8])
def test_volume_partition_conserves_and_balances_elements(ranks):
    result, rows = run_benchmark(ranks)

    assert result.returncode == 0, result.stdout
    assert len(rows) == ranks, result.stdout
    assert {row["rank"] for row in rows} == set(range(ranks))
    assert {row["ranks"] for row in rows} == {ranks}
    assert {row["before"] for row in rows} == {1660}
    assert sum(row["after"] for row in rows) == 1660

    average = 1660 / ranks
    imbalance = max(row["after"] for row in rows) / average - 1
    assert imbalance <= 0.03


def test_seed_is_deterministic():
    first, first_rows = run_benchmark(4, seed=42)
    second, second_rows = run_benchmark(4, seed=42)

    assert first.returncode == second.returncode == 0
    first_counts = sorted((row["rank"], row["after"]) for row in first_rows)
    second_counts = sorted((row["rank"], row["after"]) for row in second_rows)
    assert first_counts == second_counts


def test_cut_objective():
    result, rows = run_benchmark(4, objective="cut")

    assert result.returncode == 0, result.stdout
    assert "METIS partition: objective=cut" in result.stdout
    assert sum(row["after"] for row in rows) == 1660


def test_invalid_objective_is_rejected_on_every_rank():
    result, rows = run_benchmark(2, objective="bogus")

    # OpenSees currently exits zero for Tcl script errors, so verify the
    # interpreter diagnostics rather than the process status.
    assert result.returncode == 0
    assert not rows
    assert result.stdout.count("partition objective must be cut or volume") == 2
