"""End-to-end CLI, Python ABI, metadata matching and output schema checks."""
import copy
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


def tsv(text):
  rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
  assert all(None not in row and all(v is not None for v in row.values()) for row in rows), text
  return rows


values = {"A": 1.2, "B": 2.8, "C": 2.1, "D": 4.5}
predictor = {"A": 0.5, "B": 1.5, "C": 0.9, "D": 2.0}
classes = {"A": "low", "B": "high", "C": "low", "D": "high"}
tree_file, data_file = source / "examples/tree.nwk", source / "examples/traits.tsv"

with Tree(tree_file.read_text()) as tree:
  bm = tree.brownian(values)
  disc = tree.discrete(classes, states=["high", "low"])
  pars = tree.parsimony(classes)
  gls = tree.pgls(values, [predictor])
  assert len(tree.pic(values)) == 3
  assert pars["score"] == 2
  assert len(gls["coefficients"]) == 2
  assert tree.newick.endswith(";")
  for copier in (copy.copy, copy.deepcopy):
    try:
      copier(tree)
      raise AssertionError("Owning tree was copied")
    except TypeError:
      pass
  with tree.midpoint_root() as rooted:
    assert {n["label"] for n in rooted.nodes if n["is_tip"]} == set(values)
  extended, added = tree.extend_tips()
  with extended:
    assert added > 0
  try:
    tree.brownian({"A": 1})
    raise AssertionError("Missing metadata accepted")
  except ValueError:
    pass

try:
  tree.pic(values)
  raise AssertionError("Closed tree accepted")
except ValueError:
  pass

with Tree.upgma({"A": [0, 0], "B": [2, 0], "C": [10, 0]}) as inferred:
  assert len([n for n in inferred.nodes if n["is_tip"]]) == 3

with tempfile.TemporaryDirectory(prefix="interfaces-", dir=source / "build") as directory:
  work = Path(directory)
  for command, column in (("brownian", "value"), ("discrete", "class"), ("parsimony", "class"),
                           ("pic", "value"), ("pgls", "value")):
    prefix = work / command
    args = [command, "--tree", tree_file, "--data", data_file, "--trait", column, "--prefix", prefix]
    if command == "pgls":
      args += ["--predictor", "predictor"]
    run(*args)
    rows = tsv(prefix.with_suffix(".nodes.tsv").read_text())
    summary = {r["metric"]: r["value"] for r in tsv(prefix.with_suffix(".summary.tsv").read_text())}
    assert len(rows) == 7 and summary["method"] == command
    assert prefix.with_suffix(".tree.nwk").read_text().strip().endswith(";")
    if command == "brownian":
      assert math.isclose(float(summary["root_mean"]), bm["root_mean"], rel_tol=1e-12)
    if command == "discrete":
      assert math.isclose(float(summary["log_likelihood"]), disc["log_likelihood"], rel_tol=1e-12)
      assert all(math.isclose(float(row["p_0"]) + float(row["p_1"]), 1) for row in rows)
    if command == "parsimony":
      assert float(summary["parsimony_score"]) == pars["score"]
    if command == "pgls":
      coeff = tsv(prefix.with_suffix(".coefficients.tsv").read_text())
      assert math.isclose(float(coeff[1]["coefficient"]), gls["coefficients"][1], rel_tol=1e-12)
    if command == "pic":
      assert len(tsv(prefix.with_suffix(".contrasts.tsv").read_text())) == 3
  run("parsimony", "--tree", tree_file, "--data", data_file, "--trait", "ordered", "--ordered", "--states", 3)
  run("discrete", "--tree", tree_file, "--data", data_file, "--trait", "class", "--fit-rate")
  run("inspect", "--tree", tree_file)
  prepared = run("prepare", "--tree", tree_file, "--midpoint", "--extend-tips")
  with Tree(prepared.stdout) as transformed:
    assert len([n for n in transformed.nodes if n["is_tip"]]) == 4
  run("upgma", "--data", data_file, "--trait", "value", "--predictor", "predictor")
  lengthless, constant = work / "lengthless.nwk", work / "constant.tsv"
  lengthless.write_text("(A,B);")
  constant.write_text("taxon\tclass\nB\tred\nA\tred\n")
  valid = run("parsimony", "--tree", lengthless, "--data", constant, "--trait", "class")
  assert "parsimony_score\t0" in valid.stderr
  assert all(row["state"] == "red" for row in tsv(valid.stdout))
  run("brownian", "--tree", lengthless, "--data", constant, "--trait", "class", success=False)
  bad = work / "bad.tsv"
  bad.write_text("taxon\tvalue\nA\t1\nA\t2\nC\t3\nD\t4\n")
  run("brownian", "--tree", tree_file, "--data", bad, success=False)
  bad.write_text("taxon\tvalue\nA\t1\nB\t2\nC\t3\nE\t4\n")
  run("brownian", "--tree", tree_file, "--data", bad, success=False)
  bad.write_text("taxon\tclass\nA\tlow\nB\tNA\nC\tlow\nD\thigh\n")
  run("discrete", "--tree", tree_file, "--data", bad, "--trait", "class")
  run("discrete", "--tree", tree_file, "--data", data_file, "--trait", "ordered", "--states", 2, success=False)
  run("unknown-command", success=False)
  run("--help")
print("CLI and Python interface checks passed")
