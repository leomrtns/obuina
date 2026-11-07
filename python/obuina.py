"""Low-level ctypes binding to libobuina; no third-party Python packages.

Set OBUINA_LIBRARY to the shared library's full path when installed outside
the project. Tree owns its C pointer; use a context manager or call close(). All
returned results are Python copies, so they survive destruction of the tree.
SPDX-License-Identifier: GPL-3.0-or-later
"""
from __future__ import annotations
import ctypes as C
import ctypes.util
import os
import math
from pathlib import Path

_size = C.c_size_t
_id = C.c_uint32
_doublep = C.POINTER(C.c_double)


class _Brownian(C.Structure):
  _fields_ = [("nodes", _size), ("tips", _size), ("root_mean", C.c_double),
              ("sigma2", C.c_double), ("log_likelihood", C.c_double),
              ("mean", _doublep), ("variance", _doublep)]


class _BrownianBatch(C.Structure):
  _fields_ = [("traits", _size), ("results", C.POINTER(C.POINTER(_Brownian)))]


class _Discrete(C.Structure):
  _fields_ = [("nodes", _size), ("states", _size), ("rate", C.c_double),
              ("log_likelihood", C.c_double), ("probabilities", _doublep)]


class _Parsimony(C.Structure):
  _fields_ = [("nodes", _size), ("states", _size), ("score", C.c_double),
              ("assignment", C.POINTER(_size)), ("optimal", C.POINTER(C.c_ubyte))]


class _PIC(C.Structure):
  _fields_ = [("count", _size), ("node", C.POINTER(_id)),
              ("contrast", _doublep), ("variance", _doublep)]


class _PGLS(C.Structure):
  _fields_ = [("tips", _size), ("parameters", _size), ("coefficients", _doublep),
              ("standard_errors", _doublep), ("sigma2", C.c_double), ("log_likelihood", C.c_double)]


_lib = None


def _library():
  global _lib
  if _lib is not None:
    return _lib
  path = os.environ.get("OBUINA_LIBRARY")
  if not path:
    base = Path(__file__).resolve().parents[1]
    for candidate in (base / "build/libobuina.so", base / "build/libobuina.dylib"):
      if candidate.is_file():
        path = str(candidate)
        break
  path = path or ctypes.util.find_library("obuina")
  if not path:
    raise OSError("Build obuina or set OBUINA_LIBRARY to its shared library")
  lib = C.CDLL(path)
  signatures = {
    "bt_version": ([], C.c_char_p),
    "crpx_status_string": ([C.c_int], C.c_char_p),
    "crpx_tree_from_newick": ([C.c_char_p, C.POINTER(C.c_void_p), C.POINTER(_size)], C.c_int),
    "crpx_tree_free": ([C.c_void_p], None),
    "crpx_tree_node_count": ([C.c_void_p], _size),
    "crpx_tree_first_child": ([C.c_void_p, _id], _id),
    "crpx_tree_label": ([C.c_void_p, _id], C.c_char_p),
    "crpx_tree_to_newick": ([C.c_void_p, C.POINTER(C.c_void_p)], C.c_int),
    "crpx_tree_string_free": ([C.c_void_p], None),
    "bt_tree_midpoint_root": ([C.c_void_p, C.POINTER(C.c_void_p)], C.c_int),
    "bt_tree_extend_tips": ([C.c_void_p, C.POINTER(C.c_void_p), _doublep], C.c_int),
    "bt_tree_upgma": ([C.POINTER(C.c_char_p), _size, _size, _doublep, C.POINTER(C.c_void_p)], C.c_int),
    "bt_brownian_fit": ([C.c_void_p, _doublep, C.POINTER(C.POINTER(_Brownian))], C.c_int),
    "bt_brownian_free": ([C.POINTER(_Brownian)], None),
    "bt_brownian_fit_many": ([C.c_void_p, _doublep, _size, C.POINTER(C.POINTER(_BrownianBatch))], C.c_int),
    "bt_brownian_batch_free": ([C.POINTER(_BrownianBatch)], None),
    "bt_discrete_asr": ([C.c_void_p, C.POINTER(_size), _size, C.c_double, C.c_bool,
                         C.POINTER(C.POINTER(_Discrete))], C.c_int),
    "bt_discrete_free": ([C.POINTER(_Discrete)], None),
    "bt_discrete_asr_likelihoods": ([C.c_void_p, _doublep, _size, C.c_double, C.c_bool,
                                     C.POINTER(C.POINTER(_Discrete))], C.c_int),
    "bt_parsimony_asr": ([C.c_void_p, C.POINTER(_size), _size, C.c_bool,
                          C.POINTER(C.POINTER(_Parsimony))], C.c_int),
    "bt_parsimony_free": ([C.POINTER(_Parsimony)], None),
    "bt_pic": ([C.c_void_p, _doublep, C.POINTER(C.POINTER(_PIC))], C.c_int),
    "bt_pic_free": ([C.POINTER(_PIC)], None),
    "bt_pgls": ([C.c_void_p, _doublep, _doublep, _size, C.POINTER(C.POINTER(_PGLS))], C.c_int),
    "bt_pgls_free": ([C.POINTER(_PGLS)], None),
  }
  for name, (args, result) in signatures.items():
    function = getattr(lib, name)
    function.argtypes, function.restype = args, result
  _lib = lib
  return lib


def _check(status):
  if status:
    message = _library().crpx_status_string(status).decode()
    raise ValueError(f"obuina: {message}; see docs/methods.md for restrictions")


class Tree:
  def __init__(self, newick: str):
    if "\0" in newick:
      raise ValueError("Newick contains a NUL byte")
    self._pointer = C.c_void_p()
    offset = _size()
    status = _library().crpx_tree_from_newick(newick.encode(), C.byref(self._pointer), C.byref(offset))
    if status:
      raise ValueError(f"Newick parse failed at byte {offset.value}")

  @classmethod
  def _adopt(cls, pointer):
    tree = cls.__new__(cls)
    tree._pointer = pointer
    return tree

  def _alive(self):
    if not self._pointer:
      raise ValueError("Tree has been closed")

  def close(self):
    if getattr(self, "_pointer", None):
      _library().crpx_tree_free(self._pointer)
      self._pointer = C.c_void_p()

  def __enter__(self):
    self._alive()
    return self

  def __copy__(self):
    raise TypeError("Tree owns C memory; create a new Tree(tree.newick) to copy it")

  def __deepcopy__(self, memo):
    return self.__copy__()

  def __exit__(self, *_):
    self.close()

  def __del__(self):
    self.close()

  @property
  def nodes(self):
    self._alive()
    lib = _library()
    result = []
    for i in range(lib.crpx_tree_node_count(self._pointer)):
      label = lib.crpx_tree_label(self._pointer, i)
      result.append({"node_id": i, "label": label.decode() if label else None,
                     "is_tip": lib.crpx_tree_first_child(self._pointer, i) == 2**32 - 1})
    return result

  @property
  def newick(self):
    self._alive()
    pointer = C.c_void_p()
    _check(_library().crpx_tree_to_newick(self._pointer, C.byref(pointer)))
    try:
      return C.string_at(pointer).decode()
    finally:
      _library().crpx_tree_string_free(pointer)

  def _observations(self, values, kind=C.c_double, transform=float, allow_missing=False):
    nodes = self.nodes
    names = [node["label"] for node in nodes if node["is_tip"]]
    if any(name is None for name in names) or len(set(names)) != len(names):
      raise ValueError("Tree requires unique named tips")
    if set(values) - set(names):
      raise ValueError("Metadata contains taxa absent from the tree")
    if not allow_missing and set(values) != set(names):
      raise ValueError("Metadata and tree require the same taxon names; use allow_missing=True for omitted tips")
    array = (kind * len(nodes))()
    for node in nodes:
      if node["is_tip"]:
        array[node["node_id"]] = transform(values.get(node["label"]))
    return array

  def _continuous(self, values, allow_missing):
    def convert(value):
      if value is None:
        if allow_missing:
          return math.nan
        raise ValueError("Unknown continuous observations require allow_missing=True")
      number = float(value)
      if math.isnan(number) and not allow_missing:
        raise ValueError("NaN observations require allow_missing=True")
      if math.isinf(number):
        raise ValueError("Infinite observations are invalid")
      return number
    return self._observations(values, transform=convert, allow_missing=allow_missing)

  def midpoint_root(self):
    """Return a new tree; node IDs may change."""
    self._alive()
    pointer = C.c_void_p()
    _check(_library().bt_tree_midpoint_root(self._pointer, C.byref(pointer)))
    return self._adopt(pointer)

  def extend_tips(self):
    """Return (new tree, added length), with internal edges fixed."""
    self._alive()
    pointer, added = C.c_void_p(), C.c_double()
    _check(_library().bt_tree_extend_tips(self._pointer, C.byref(pointer), C.byref(added)))
    return self._adopt(pointer), added.value

  @classmethod
  def upgma(cls, metadata):
    """Build a descriptive tree from taxon -> sequence of scaled numeric features."""
    names = list(metadata)
    if len(names) < 2 or any(not isinstance(name, str) or not name or "\0" in name for name in names):
      raise ValueError("UPGMA requires >=2 nonempty string taxon names")
    features = len(metadata[names[0]])
    if not features or any(len(metadata[name]) != features for name in names):
      raise ValueError("UPGMA feature rows must have the same positive length")
    labels = (C.c_char_p * len(names))(*(name.encode() for name in names))
    values = (C.c_double * (len(names) * features))(*(float(x) for name in names for x in metadata[name]))
    pointer = C.c_void_p()
    _check(_library().bt_tree_upgma(labels, len(names), features, values, C.byref(pointer)))
    return cls._adopt(pointer)

  def _brownian_dict(self, fit, values):
    nodes = self.nodes
    for node in nodes:
      i = node["node_id"]
      observed = node["is_tip"] and values.get(node["label"]) is not None
      if observed:
        observed = not math.isnan(float(values[node["label"]]))
      node.update(mean=fit.mean[i], variance=fit.variance[i], observed=observed)
    return {"root_mean": fit.root_mean, "sigma2": fit.sigma2, "log_likelihood": fit.log_likelihood,
            "observed_tips": fit.tips, "nodes": nodes}

  def brownian(self, values, allow_missing=False):
    """Fit observed-data likelihood; predict missing tips with conditional variance."""
    array, result = self._continuous(values, allow_missing), C.POINTER(_Brownian)()
    _check(_library().bt_brownian_fit(self._pointer, array, C.byref(result)))
    try:
      return self._brownian_dict(result.contents, values)
    finally:
      _library().bt_brownian_free(result)

  def brownian_many(self, traits, allow_missing=False):
    """Independent fits: mapping trait name -> taxon/value mapping; no cross-trait borrowing."""
    names = list(traits)
    if not names:
      raise ValueError("At least one trait is required")
    columns = [self._continuous(traits[name], allow_missing) for name in names]
    nodes = len(self.nodes)
    values = (C.c_double * (nodes * len(names)))(
      *(columns[j][i] for i in range(nodes) for j in range(len(names))))
    result = C.POINTER(_BrownianBatch)()
    _check(_library().bt_brownian_fit_many(self._pointer, values, len(names), C.byref(result)))
    try:
      return {name: self._brownian_dict(result.contents.results[j].contents, traits[name])
              for j, name in enumerate(names)}
    finally:
      _library().bt_brownian_batch_free(result)

  def _states(self, values, states, allow_missing=False):
    states = list(states) if states is not None else sorted({x for x in values.values() if x is not None}, key=str)
    if not states or None in states or len(set(states)) != len(states):
      raise ValueError("Supply a nonempty unique state alphabet; None is reserved for unknown")
    codes = {state: i for i, state in enumerate(states)}
    missing = _size(-1).value
    array = self._observations(values, _size, lambda x: missing if x is None else codes[x], allow_missing)
    return states, array

  def discrete(self, values, states=None, rate=1.0, fit_rate=False, allow_missing=False):
    states, array = self._states(values, states, allow_missing)
    result = C.POINTER(_Discrete)()
    _check(_library().bt_discrete_asr(self._pointer, array, len(states), rate, fit_rate, C.byref(result)))
    try:
      return self._discrete_dict(result.contents, states, values)
    finally:
      _library().bt_discrete_free(result)

  def _discrete_dict(self, fit, states, values):
    nodes = self.nodes
    for node in nodes:
      i = node["node_id"]
      node["probabilities"] = [fit.probabilities[i * len(states) + s] for s in range(len(states))]
      node["observed"] = node["is_tip"] and values.get(node["label"]) is not None
    return {"states": states, "rate": fit.rate, "log_likelihood": fit.log_likelihood, "nodes": nodes}

  def discrete_likelihoods(self, values, states, rate=1.0, fit_rate=False, allow_missing=False):
    """Taxon -> likelihood vector p(observation|state); None means all ones.

    Weights need not sum to one. State posterior vectors should be converted using
    the prior used to generate them; see docs/missing-data-and-uncertainty.md.
    """
    states = list(states)
    if len(states) < 2 or None in states or len(set(states)) != len(states):
      raise ValueError("An explicit alphabet of >=2 unique states is required")
    self._observations(values, transform=lambda _: 0, allow_missing=allow_missing)
    nodes = self.nodes
    array = (C.c_double * (len(nodes) * len(states)))()
    for node in nodes:
      if not node["is_tip"]:
        continue
      vector = values.get(node["label"])
      vector = [1.0] * len(states) if vector is None else list(vector)
      if len(vector) != len(states):
        raise ValueError("Likelihood vectors must match the state alphabet")
      for j, value in enumerate(vector):
        array[node["node_id"] * len(states) + j] = float(value)
    result = C.POINTER(_Discrete)()
    _check(_library().bt_discrete_asr_likelihoods(self._pointer, array, len(states), rate, fit_rate,
                                               C.byref(result)))
    try:
      return self._discrete_dict(result.contents, states, values)
    finally:
      _library().bt_discrete_free(result)

  def parsimony(self, values, states=None, ordered=False, allow_missing=False):
    if ordered and states is None:
      raise ValueError("Ordered parsimony requires an explicit ordered state alphabet")
    states, array = self._states(values, states, allow_missing)
    result = C.POINTER(_Parsimony)()
    _check(_library().bt_parsimony_asr(self._pointer, array, len(states), ordered, C.byref(result)))
    try:
      fit, nodes = result.contents, self.nodes
      for node in nodes:
        i = node["node_id"]
        node["assignment"] = states[fit.assignment[i]]
        node["observed"] = node["is_tip"] and values.get(node["label"]) is not None
        node["optimal_states"] = [states[s] for s in range(len(states)) if fit.optimal[i * len(states) + s]]
      return {"score": fit.score, "states": states, "nodes": nodes}
    finally:
      _library().bt_parsimony_free(result)

  def pic(self, values, allow_missing=False):
    """Contrasts from observed tips; never substitute predicted states."""
    array, result = self._continuous(values, allow_missing), C.POINTER(_PIC)()
    _check(_library().bt_pic(self._pointer, array, C.byref(result)))
    try:
      fit = result.contents
      return [{"node_id": fit.node[i], "contrast": fit.contrast[i], "variance": fit.variance[i]}
              for i in range(fit.count)]
    finally:
      _library().bt_pic_free(result)

  def pgls(self, response, predictors=(), allow_missing=False):
    """Predictors are taxon->number mappings; intercept added, incomplete rows omitted."""
    y = self._continuous(response, allow_missing)
    columns = [self._continuous(values, allow_missing) for values in predictors]
    nodes = len(self.nodes)
    x = (C.c_double * (nodes * len(columns)))(*(columns[j][i] for i in range(nodes) for j in range(len(columns))))
    result = C.POINTER(_PGLS)()
    _check(_library().bt_pgls(self._pointer, y, x, len(columns), C.byref(result)))
    try:
      fit = result.contents
      return {"coefficients": list(fit.coefficients[:fit.parameters]),
              "standard_errors": list(fit.standard_errors[:fit.parameters]),
              "sigma2": fit.sigma2, "log_likelihood": fit.log_likelihood, "used_tips": fit.tips}
    finally:
      _library().bt_pgls_free(result)
