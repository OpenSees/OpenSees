import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXE = ROOT / "build" / "mpich-release" / "OpenSeesMP"
BENCHMARK = ROOT / "EXAMPLES" / "ParallelModelMP" / "metisPartitionBenchmark.tcl"
STATE_MODEL = ROOT / "EXAMPLES" / "ParallelModelMP" / "metisPartitionState.tcl"
HOMEBREW_TCL_LIBRARY = Path("/opt/homebrew/opt/tcl-tk@8/lib/tcl8.6")
ROW = re.compile(
    r"METIS_BENCHMARK rank=(?P<rank>\d+) ranks=(?P<ranks>\d+) "
    r"before=(?P<before>\d+)\s+after=(?P<after>\d+)"
)
STATE_ROW = re.compile(r"^METIS_STATE (?P<fields>.+)$", re.MULTILINE)


def executable():
    path = Path(os.environ.get("OPENSEESMP", DEFAULT_EXE))
    if not path.is_file():
        pytest.skip(f"OpenSeesMP executable not found: {path}")
    return path


def mpi_launcher():
    command = os.environ.get("MPIEXEC", "mpiexec")
    launcher = shutil.which(command)
    if launcher is None:
        pytest.skip(f"MPI launcher not found: {command}")
    return launcher


def run_model(ranks, script, *args):
    env = os.environ.copy()
    if "TCL_LIBRARY" not in env and HOMEBREW_TCL_LIBRARY.is_dir():
        env["TCL_LIBRARY"] = str(HOMEBREW_TCL_LIBRARY)
    launcher_flags = shlex.split(os.environ.get("MPIEXEC_FLAGS", ""))
    result = subprocess.run(
        [mpi_launcher(), *launcher_flags, "-n", str(ranks),
         str(executable()), str(script), *map(str, args)],
        cwd=ROOT,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=30,
    )
    return result


def run_benchmark(ranks, objective="volume", seed=42):
    result = run_model(ranks, BENCHMARK, objective, seed)
    return result, [
        {key: int(value) for key, value in match.groupdict().items()}
        for match in ROW.finditer(result.stdout)
    ]


def csv_ints(value):
    return [] if not value else [int(item) for item in value.split(",")]


def state_rows(output):
    rows = []
    for match in STATE_ROW.finditer(output):
        fields = dict(item.split("=", 1)
                      for item in match.group("fields").split())
        row = {
            "rank": int(fields["rank"]),
            "ranks": int(fields["ranks"]),
        }
        for name in ("nodes", "elements", "eleloads", "nodeloads",
                     "fixed", "constrained", "mp100", "mp102"):
            row[name] = csv_ints(fields[name])
        row["masses"] = {
            int(item.split(":", 1)[0]): float(item.split(":", 1)[1])
            for item in fields["masses"].split(",") if item
        }
        rows.append(row)
    return rows


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
    assert len(first_rows) == len(second_rows) == 4
    first_counts = sorted((row["rank"], row["after"]) for row in first_rows)
    second_counts = sorted((row["rank"], row["after"]) for row in second_rows)
    assert first_counts == second_counts


def test_cut_objective():
    result, rows = run_benchmark(4, objective="cut")

    assert result.returncode == 0, result.stdout
    assert "METIS partition: objective=cut" in result.stdout
    assert len(rows) == 4
    assert sum(row["after"] for row in rows) == 1660


def test_invalid_objective_is_rejected_on_every_rank():
    result, rows = run_benchmark(2, objective="bogus")

    # OpenSees currently exits zero for Tcl script errors, so verify the
    # interpreter diagnostics rather than the process status.
    assert result.returncode == 0
    assert not rows
    assert result.stdout.count("partition objective must be cut or volume") == 2


def test_partition_preserves_load_mass_and_constraint_ownership():
    result = run_model(4, STATE_MODEL)
    rows = state_rows(result.stdout)

    assert result.returncode == 0, result.stdout
    assert len(rows) == 4, result.stdout
    assert {row["rank"] for row in rows} == set(range(4))
    rows_by_rank = {row["rank"]: row for row in rows}

    # Each element and its element load belong to exactly one, identical rank.
    for row in rows:
        assert row["eleloads"] == row["elements"]
    assert sorted(tag for row in rows for tag in row["elements"]) == [1, 2, 3, 4]

    expected_nodes = {1, 2, 3, 4, 5, 100, 101, 102, 103}
    mass_owners = {}
    for tag in expected_nodes:
        owners = [row["rank"] for row in rows
                  if row["masses"].get(tag, 0.0) > 0.0]
        assert len(owners) == 1, (tag, rows, result.stdout)
        owner = owners[0]
        assert rows_by_rank[owner]["masses"][tag] == pytest.approx(float(tag))
        mass_owners[tag] = owner

    # The mesh-to-floating MP pair has the same owner and is complete on every
    # rank that retains the mesh interface node.
    assert mass_owners[3] == mass_owners[100]
    ranks_with_3 = {row["rank"] for row in rows if 3 in row["nodes"]}
    ranks_with_100 = {row["rank"] for row in rows if 100 in row["nodes"]}
    assert ranks_with_100 == ranks_with_3
    for row in rows:
        assert row["mp100"] == ([3] if row["rank"] in ranks_with_3 else [])

    # A floating-to-floating MP pair and a fixed floating node each survive on
    # exactly one rank, with their constraint and applied nodal load intact.
    assert mass_owners[101] == mass_owners[102]
    for row in rows:
        has_floating_pair = row["rank"] == mass_owners[101]
        assert (101 in row["nodes"]) == has_floating_pair
        assert (102 in row["nodes"]) == has_floating_pair
        assert row["mp102"] == ([101] if has_floating_pair else [])

    fixed_rows = [row for row in rows if 103 in row["fixed"]]
    assert [row["rank"] for row in fixed_rows] == [mass_owners[103]]
    assert sorted(tag for row in rows for tag in row["nodeloads"]) == [100, 103]
    assert all(tag in row["nodes"]
               for row in rows for tag in row["nodeloads"])
