"""Missing observations, independent dimensions and uncertainty through CLI/Python."""
import csv
import io
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

executable, library, source = map(Path, sys.argv[1:])
os.environ["OBUINA_LIBRARY"] = str(library)
sys.path.insert(0, str(source / "python"))
from obuina import Tree


def run(*args, success=True):
  result = subprocess.run([str(executable), *map(str, args)], capture_output=True, text=True)
  assert (result.returncode == 0) == success, (args, result.stdout, result.stderr)
  return result


def rows(path):
  return list(csv.DictReader(Path(path).open(), delimiter="\t"))


def metrics(prefix):
  return {r["metric"]: r["value"] for r in rows(str(prefix) + ".summary.tsv")}


def fails(function):
  try:
    function()
  except (ValueError, KeyError):
    return
  raise AssertionError("Invalid/missing observations unexpectedly accepted")


tree_file, data_file = source / "examples/tree.nwk", source / "examples/traits.tsv"
y = {"A": 1.2, "B": 2.8, "C": 2.1, "D": 4.5}
x = {"A": 0.5, "B": 1.5, "C": 0.9, "D": 2.0}
x2 = {"A": 0, "B": 2, "C": 1, "D": 2}
with Tree(tree_file.read_text()) as tree:
  known = {"A": 1.2, "C": 2.1}
  fit = tree.brownian(known, allow_missing=True)
  explicit = tree.brownian({"A": 1.2, "B": None, "C": 2.1, "D": math.nan}, allow_missing=True)
  assert fit["observed_tips"] == 2
  assert fit["root_mean"] == explicit["root_mean"]
  for node in fit["nodes"]:
    if node["is_tip"]:
      assert node["observed"] == (node["label"] in known)
      assert (node["variance"] == 0) == node["observed"]
  fails(lambda: tree.brownian(known))
  fails(lambda: tree.brownian({"A": 1.2, "B": None, "C": 2.1, "D": 4.5}))
  fails(lambda: tree.brownian({"A": 1}, allow_missing=True))
  fails(lambda: tree.brownian({"A": 1, "C": 2, "TYPO": 4}, allow_missing=True))
  fails(lambda: tree.brownian({"A": math.inf, "C": 2}, allow_missing=True))
  many = tree.brownian_many({"first": known, "second": {"B": 5, "C": 8, "D": 2}}, allow_missing=True)
  assert many["first"]["observed_tips"] == 2 and many["second"]["observed_tips"] == 3
  assert many["first"]["log_likelihood"] == fit["log_likelihood"]
  fails(lambda: tree.brownian_many({"first": known, "empty": {}}, allow_missing=True))
  assert len(tree.pic(known, allow_missing=True)) == 1
  assert len(tree.pic({"A": 1.2, "C": 2.1, "D": 4.5}, allow_missing=True)) == 2
  missing_gls = tree.pgls({"A": 1.2, "B": 2.8, "C": 2.1}, [x], allow_missing=True)
  assert missing_gls["used_tips"] == 3
  assert tree.pgls(y, [x, x2])["used_tips"] == 4
  fails(lambda: tree.pgls(known, [x], allow_missing=True))
  classes = {"A": "low", "C": "high"}
  disc = tree.discrete(classes, states=["low", "high"], allow_missing=True)
  disc_none = tree.discrete({"A": "low", "B": None, "C": "high", "D": None}, states=["low", "high"])
  assert disc["log_likelihood"] == disc_none["log_likelihood"]
  for a, b in zip(disc["nodes"], disc_none["nodes"]):
    assert a["probabilities"] == b["probabilities"]
  assert tree.parsimony(classes, states=["low", "high"], allow_missing=True)["score"] == 1
  fails(lambda: tree.discrete(classes, states=["low", "high"]))
  weights = {"A": [0.8, 0.2], "C": [0.1, 0.9]}
  soft = tree.discrete_likelihoods(weights, ["low", "high"], allow_missing=True)
  scaled = tree.discrete_likelihoods({"A": [8, 2], "C": [0.1, 0.9]}, ["low", "high"], allow_missing=True)
  assert math.isclose(scaled["log_likelihood"] - soft["log_likelihood"], math.log(10), abs_tol=1e-12)
  for a, b in zip(soft["nodes"], scaled["nodes"]):
    assert all(math.isclose(x, y, abs_tol=1e-12) for x, y in zip(a["probabilities"], b["probabilities"]))
  for vector in ([0, 0], [-1, 2], [math.nan, 1], [math.inf, 1], [1, 2, 3]):
    fails(lambda: tree.discrete_likelihoods({"A": vector}, [0, 1], allow_missing=True))

with tempfile.TemporaryDirectory(prefix="missing-", dir=source / "build") as directory:
  work = Path(directory)
  data = work / "partial.tsv"
  data.write_text("taxon\tvalue\tsecond\tpredictor\tclass\nA\t1.2\t10\t0.5\tlow\n"
                  "B\tNA\t14\t1.5\t?\nC\t2.1\t?\t0.9\thigh\n")
  prefix = work / "missing"
  run("brownian", "--tree", tree_file, "--data", data, "--prefix", prefix, success=False)
  run("brownian", "--tree", tree_file, "--data", data, "--missing", "estimate", "--prefix", prefix)
  result = metrics(prefix)
  assert result["observed_tips"] == "2" and result["missing_tips"] == "2" and result["used_tips"] == "2"
  assert math.isclose(float(result["root_mean"]), fit["root_mean"], abs_tol=1e-12)
  tips = {r["label"]: r for r in rows(str(prefix) + ".nodes.tsv") if r["is_tip"] == "1"}
  assert tips["B"]["observed"] == "0" and tips["D"]["observed"] == "0"
  assert float(tips["D"]["variance"]) > 0 and tips["D"]["used_for_fit"] == "0"
  for command in ("pic", "discrete", "parsimony"):
    args = [command, "--tree", tree_file, "--data", data, "--missing", "estimate", "--prefix", prefix]
    if command in ("discrete", "parsimony"):
      args += ["--trait", "class"]
    run(*args)
  prefix = work / "batch"
  run("brownian", "--tree", tree_file, "--data", data, "--traits", "value,second",
      "--missing", "estimate", "--prefix", prefix)
  manifest = rows(str(prefix) + ".traits.tsv")
  assert [r["trait"] for r in manifest] == ["value", "second"]
  assert all(r["status"] == "ok" for r in manifest)
  assert metrics(manifest[0]["prefix"])["trait"] == "value"
  assert metrics(manifest[1]["prefix"])["trait"] == "second"
  run("brownian", "--tree", tree_file, "--data", data, "--traits", "value,second", success=False)
  run("brownian", "--tree", tree_file, "--data", data, "--traits", "value,value", "--prefix", prefix,
      success=False)
  run("brownian", "--tree", tree_file, "--data", data, "--traits", "value,second", "--trait", "value",
      "--prefix", prefix, success=False)
  run("brownian", "--tree", tree_file, "--data", data, "--traits", "absent,value", "--missing", "estimate",
      "--prefix", prefix, success=False)
  assert [r["status"] for r in rows(str(prefix) + ".traits.tsv")] == ["error", "ok"]
  run("discrete", "--tree", tree_file, "--data", data_file, "--traits", "class,ordered", "--prefix", prefix)
  run("pgls", "--tree", tree_file, "--data", data_file, "--predictors", "predictor,ordered", "--prefix", prefix)
  coefficients = rows(str(prefix) + ".coefficients.tsv")
  assert [r["parameter"] for r in coefficients] == ["intercept", "predictor", "ordered"]
  with Tree(tree_file.read_text()) as tree:
    expected = tree.pgls(y, [x, x2])
  for row, beta in zip(coefficients, expected["coefficients"]):
    assert math.isclose(float(row["coefficient"]), beta, abs_tol=1e-11)
  data.write_text("taxon\tvalue\tpredictor\nA\t1.2\t0.5\nB\t2.8\t1.5\nC\t2.1\tNA\nD\t4.5\t2.0\n")
  run("pgls", "--tree", tree_file, "--data", data, "--predictor", "predictor", "--missing", "estimate",
      "--prefix", prefix)
  assert metrics(prefix)["observed_tips"] == "4" and metrics(prefix)["used_tips"] == "3"
  tips = {r["label"]: r for r in rows(str(prefix) + ".nodes.tsv") if r["is_tip"] == "1"}
  assert tips["C"]["observed"] == "1" and tips["C"]["used_for_fit"] == "0"
  run("upgma", "--data", data_file, "--predictors", "predictor,ordered")
  run("upgma", "--data", data, "--predictors", "predictor", "--missing", "estimate", success=False)
  for text in ("taxon\tvalue\nA\t1\nC\t2\nTYPO\t3\n", "taxon\tvalue\nA\t1\nA\t3\nC\t2\n"):
    data.write_text(text)
    run("brownian", "--tree", tree_file, "--data", data, "--missing", "estimate", success=False)
print("Missing, multi-column and uncertain-observation interface checks passed")
