"""Small reproducible scaling benchmark; run from the project after a Release build.

PYTHONPATH=python python3 tests/benchmark.py --tips 1000 10000
Times include binding calls, result copying and per-call allocation; tree parsing
is outside timing. This is descriptive and imposes no machine-specific threshold.
"""
import argparse
import math
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from obuina import Tree


def balanced(begin, end):
  if end - begin == 1:
    return f"T{begin}:1"
  middle = (begin + end) // 2
  return f"({balanced(begin, middle)},{balanced(middle, end)}):1"


parser = argparse.ArgumentParser()
parser.add_argument("--tips", nargs="+", type=int, default=[1000, 10000])
args = parser.parse_args()
print("tips\tmethod\tseconds")
for tips in args.tips:
  if tips < 2:
    parser.error("Need at least two tips")
  continuous = {f"T{i}": math.sin(i * 0.7) + (i % 17) / 3 for i in range(tips)}
  discrete = {f"T{i}": i % 3 for i in range(tips)}
  with Tree(balanced(0, tips) + ";") as tree:
    for name, calculation in (("brownian", lambda: tree.brownian(continuous)),
                               ("pic", lambda: tree.pic(continuous)),
                               ("discrete_fixed", lambda: tree.discrete(discrete, states=[0, 1, 2])),
                               ("parsimony", lambda: tree.parsimony(discrete, states=[0, 1, 2]))):
      start = time.perf_counter()
      result = calculation()
      elapsed = time.perf_counter() - start
      del result
      print(f"{tips}\t{name}\t{elapsed:.6f}")
