"""Run from project root: PYTHONPATH=python python3 examples/missing_and_uncertain.py.

Missing continuous observations and omitted taxa require allow_missing=True.
Discrete vectors encode p(observation | latent state), not classifier posteriors.
SPDX-License-Identifier: GPL-3.0-or-later
"""
from obuina import Tree

newick = "((A:1,B:2)AB:1,(C:0.8,D:1.4)CD:1.5)root;"
with Tree(newick) as tree:
  # B is explicitly unknown and D is omitted. Fit only A and C, predict all tips.
  values = {"A": 1.2, "B": None, "C": 2.1}
  fit = tree.brownian(values, allow_missing=True)
  print("Brownian observed tips:", fit["observed_tips"])
  for node in fit["nodes"]:
    if node["is_tip"]:
      print(node["label"], "mean:", node["mean"], "variance:", node["variance"])

  columns = {"value": values, "second": {"A": 10.0, "B": 14.0, "C": None}}
  fits = tree.brownian_many(columns, allow_missing=True)
  print("Independent column fits:", list(fits))

  # Alphabet order fixes the interpretation of every vector component.
  weights = {"A": [0.9, 0.1], "B": [0.2, 0.8], "C": [1.0, 1.0]}
  uncertain = tree.discrete_likelihoods(weights, states=["low", "high"], allow_missing=True)
  print("Discrete observation likelihood:", uncertain["log_likelihood"])
  for node in uncertain["nodes"]:
    if node["is_tip"]:
      print(node["label"], "state probabilities:", node["probabilities"])
