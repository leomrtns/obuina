"""Run after building: PYTHONPATH=python python3 examples/python_example.py"""
from obuina import Tree

values = {"A": 1.2, "B": 2.8, "C": 2.1, "D": 4.5}
predictor = {"A": 0.5, "B": 1.5, "C": 0.9, "D": 2.0}
classes = {"A": "low", "B": "high", "C": "low", "D": "high"}

with Tree("((A:1,B:2):1,(C:0.8,D:1.4):1.5);") as tree:
  print("Brownian fit:", tree.brownian(values))
  print("Discrete ASR:", tree.discrete(classes, states=["low", "high"]))
  print("Parsimony:", tree.parsimony(classes))
  print("Independent contrasts:", tree.pic(values))
  print("PGLS:", tree.pgls(values, [predictor]))
