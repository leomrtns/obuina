# Validation record

The first development checks ran on Linux with GCC 13.3.0 and Python 3.12.8 on 7 November 2026. The curupixa checkout was `268040a6f359b7b3d98f346211076559d667c6b6`. The dependency source was left unchanged.

## Initial 0.1.0 checks

The Release CTest suite passed **5/5** test groups:

| Group | Independent checks and covered cases |
|---|---|
| `trees` | Every tip-pair distance before/after midpoint rooting; expected midpoint heights; exact-vertex and interior-edge roots; unary stems, zero edges, label preservation, UPGMA and restricted terminal extension |
| `continuous` | Brownian profile likelihood against a separately constructed dense tip covariance; all internal ASR means/variances against inversion of a dense Gaussian precision matrix; PIC; closed-form and multiple-predictor GLS; invalid and degenerate inputs |
| `discrete` | Enumeration of all joint assignments to check likelihood and full marginals; enumeration of ordered/unordered parsimony scores, global optimal sets and coherent traceback; missing observations, unary nodes, polytomies, zero edges, bounded rate fit and severe underflow |
| `interfaces` | Every CLI method, metadata/name matching, output tables, text alphabets, ordered codes, branchless and one-class parsimony, preparation, Python ABI and result parity, closed-tree and ownership-copy rejection |
| `c_example` | Public C API use, owned results and cleanup |

The discrete enumeration uses curupixa's transition-matrix interface independently of its pruning calculation. An additional 600-tip test with $10^{-200}$-length branches checks the analytically expected binary log likelihood and normalised marginals. The continuous references use explicit dense covariance/precision calculations independent of the pruning implementation.

AddressSanitizer and UndefinedBehaviourSanitizer passed **4/4 native CTest groups** in a Debug build. The Python/CLI interface harness also passed when using the sanitised CLI with the Release shared library for the Python binding. The Python interpreter itself was not rebuilt with sanitisation. LeakSanitizer cannot run under this environment's ptrace restrictions, so these checks used `ASAN_OPTIONS=detect_leaks=0`; they do not establish absence of leaks.

A project-local installation in `stage/` passed an installed CLI invocation and a separate CMake consumer using `find_package(obuina CONFIG REQUIRED)`, installed headers and the exported target. The Python example ran successfully. No external R/PastML cross-package comparison has yet been run; [VALIDATE-01](limitations.md) tracks matching model conventions before those comparisons.

## Reproduce

From the project directory:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
ctest --test-dir build --output-on-failure
PYTHONPATH=python python3 examples/python_example.py

cmake -S . -B build-sanitise -DCMAKE_BUILD_TYPE=Debug -DBT_SANITISE=ON
cmake --build build-sanitise -j 4
ASAN_OPTIONS=detect_leaks=0 ctest --test-dir build-sanitise --output-on-failure
ASAN_OPTIONS=detect_leaks=0 python3 tests/test_interfaces.py "$PWD/build-sanitise/obuina" "$PWD/build/libobuina.so" "$PWD"
```

The normal build tests the Python binding. The sanitised CTest configuration skips that binding test because loading an ASan library into an ordinary Python interpreter requires additional runtime configuration; the separate invocation above exercises the sanitised CLI instead.

## Descriptive scaling check

The [benchmark](../tests/benchmark.py) generates balanced binary trees with unit edges and deterministic continuous/three-state observations. Timings include Python input mapping, per-call allocation and copying result objects back to Python; tree generation/parsing is outside the timed section. These are single-run wall times on the development environment, with shared-host variation and no performance guarantee. They exclude rate fitting, midpoint rooting, UPGMA and dense PGLS.

| Tips | Brownian ASR | PIC | Fixed-rate three-state ASR | Three-state parsimony |
|---:|---:|---:|---:|---:|
| 1,000 | 0.0069 s | 0.0029 s | 0.0089 s | 0.0068 s |
| 10,000 | 0.1015 s | 0.0467 s | 0.1039 s | 0.0637 s |
| 100,000 | 0.5806 s | 0.3216 s | 0.9121 s | 0.6551 s |

```sh
PYTHONPATH=python python3 tests/benchmark.py --tips 1000 10000 100000
```

These timings and complexity exceptions describe 0.1.0. Its wide-polytomy midpoint rebuilding could be quadratic. Version 0.2.0 adopts the upstream builder to remove that append bottleneck; dense PGLS still requires quadratic storage and cubic solve time. Current costs appear in [limitations.md](limitations.md) and [curupixa-issues.md](curupixa-issues.md).

## 0.2.0 adoption and missing-data checks

The earlier results and timings describe the initial 0.1.0 implementation. Version 0.2.0 adopts curupixa `a1adfc4258fe9c6bc771464934db9d2ef6d5d397`; its tail-aware builder removes the former broad-polytomy append bottleneck, and its prepared likelihood/reconstruction/continuous APIs replace duplicated kernels. The historical timing table has not been relabelled as a new benchmark.

The updated Release suite passed **6/6 groups**: `trees`, `continuous`, `discrete`, `interfaces`, `missing_interfaces` and `c_example`. It covers missing-tip likelihoods against the original-root observed covariance, all-node predictions against dense conditional Gaussian calculations, PIC contraction and node correspondence, complete-case/multiple-predictor GLS, uncertain discrete likelihoods against exhaustive enumeration, scaling/log-likelihood offsets, independent batch ownership and distinct missingness masks, and CLI/Python policy/provenance. The 4,096-tip star regression verifies the new builder path without imposing a machine-dependent timing threshold.

The Debug ASan/UBSan suite passed **4/4 native groups** with `ASAN_OPTIONS=detect_leaks=0`. Both interface harnesses also passed against the sanitised CLI and Release Python library; LeakSanitizer remains unavailable in this environment. Strict C17 warning checks passed for all application/library sources. The updated project-local installation, independent installed CMake consumer and installed CLI version check passed at 0.2.0. Curupixa's source remained unchanged.

The [missing-data guide](missing-data-and-uncertainty.md) provides reproducible examples; [examples/missing_and_uncertain.py](../examples/missing_and_uncertain.py) exercises the new binding methods. To run the additional interface harness against the sanitised CLI:

```sh
ASAN_OPTIONS=detect_leaks=0 python3 tests/test_missing_interfaces.py "$PWD/build-sanitise/obuina" "$PWD/build/libobuina.so" "$PWD"
```
