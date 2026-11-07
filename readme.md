# Obuina

The name comes from a Brazilian tree. The public C prefix `bt_` continues to mean **branching traits**.

Obuina is a C17 library and command-line tool for ancestral state reconstruction (ASR) and phylogenetic comparative methods (PCM). It models numerical measurements, categorical classes and ordered states along a supplied tree. Taxa can represent species, samples or other entities with a defensible shared branching history; sequence data are optional.

Version **0.2.0** uses the local [curupixa](../curupixa/) library for trees, Newick parsing, UPGMA, prepared likelihood evaluation, discrete/parsimony reconstruction, Brownian fitting/reconstruction and PIC. Obuina supplies the application interfaces, missing-data handling, independent trait batches and Brownian PGLS.

**Start with the [ASR and PCM tutorial](docs/asr-pcm-tutorial.md), [missing-data and uncertainty guide](docs/missing-data-and-uncertainty.md), [method assumptions](docs/methods.md), [current limitations and work queue](docs/limitations.md), and [curupixa adoption report](docs/curupixa-issues.md).**

## Implemented methods

| Command | Purpose | Trait type | Tree requirements |
|---|---|---|---|
| `brownian` | Brownian ML fit and all-node ASR means/variances | Continuous | Rooted, positive lengths; polytomies accepted |
| `discrete` | Equal-rates CTMC fit and all-node marginal probabilities | Unordered categorical | Rooted, nonnegative lengths; unknown tips accepted |
| `parsimony` | One optimal ASR and all globally optimal state sets | Unordered or ordered | Rooted representation; lengths optional; hard polytomies |
| `pic` | Standardised phylogenetic independent contrasts | Continuous | Rooted, strictly bifurcating, positive lengths |
| `pgls` | Brownian phylogenetic generalised least squares | Continuous response/predictors | Rooted, positive lengths; polytomies accepted |
| `inspect` | Validate tree and report ultrametricity | — | Complete nonnegative lengths, unique named tips |
| `prepare` | Optional midpoint rooting and terminal extension | — | Complete nonnegative lengths |
| `upgma` | Descriptive tree from Euclidean metadata distances | Complete numeric features | No input tree required |

All fitted methods accept non-ultrametric trees. Branch lengths express the amount of evolution or elapsed time assumed by the chosen model. The root supplied in Newick is used as given; choose `--midpoint` explicitly when appropriate.

## Build and test

From this directory, with a C compiler, CMake ≥3.16 and the sibling `../curupixa` source available:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
ctest --test-dir build --output-on-failure
build/obuina --help
```

Use `-DBT_CURUPIXA_SOURCE_DIR=/path/to/curupixa` for another source location. The checkout must provide the continuous-trait, reconstruction, prepared-likelihood and tree-builder APIs; the [adoption report](docs/curupixa-issues.md) records the tested commit. No external numerical or Python packages are required. Python 3 enables the interface test and example. The build produces `build/libobuina.so` and `build/obuina` on Linux. CMake embeds the required curupixa modules in the shared library; it does not change curupixa's source or its existing build.

For a project-local installation:

```sh
cmake --install build --prefix "$PWD/stage"
stage/bin/obuina --version
```

The package installs CMake configuration, public C headers and the Python binding. Linux is the tested platform; portability work is tracked in the limitations document.

## Metadata and examples

Input is UTF-8 tab-separated text with a header. The `taxon` column matches the tree's unique tip labels exactly; row order is arbitrary. Duplicate and extra taxa are errors. Missing rows are errors by default and can be allowed explicitly with `--missing estimate`. Select one trait with `--trait`, or several independent analyses with `--traits`. Example data are in [examples/traits.tsv](examples/traits.tsv), [examples/traits-missing.tsv](examples/traits-missing.tsv) and [examples/tree.nwk](examples/tree.nwk).

```sh
mkdir -p results
build/obuina brownian --tree examples/tree.nwk --data examples/traits.tsv --prefix results/brownian
build/obuina discrete --tree examples/tree.nwk --data examples/traits.tsv --trait class --fit-rate --prefix results/classes
build/obuina parsimony --tree examples/tree.nwk --data examples/traits.tsv --trait ordered --ordered --states 3 --prefix results/ordered
build/obuina pic --tree examples/tree.nwk --data examples/traits.tsv --prefix results/pic
build/obuina pgls --tree examples/tree.nwk --data examples/traits.tsv --predictor predictor --prefix results/regression
```

Unordered classes can use text labels; their state alphabet is inferred in lexicographic order. To retain unobserved states, supply integer codes and `--states K`. Ordered parsimony requires `--states K` and codes `0..K-1`: each step costs one, including steps through unobserved intermediate states. For discrete traits, `?`, `NA` and an empty field mean unknown. Continuous markers require `--missing estimate`; supplied numeric observations must be finite.

## Missing tips and several metadata columns

```sh
build/obuina brownian --tree examples/tree.nwk --data examples/traits-missing.tsv --trait value --missing estimate --prefix results/missing-value
build/obuina brownian --tree examples/tree.nwk --data examples/traits-missing.tsv --traits value,second --missing estimate --prefix results/columns
build/obuina pgls --tree examples/tree.nwk --data examples/traits.tsv --predictors predictor,ordered --prefix results/multiple-predictors
```

Brownian fitting integrates missing tip states out of the likelihood and returns conditional predictions and variances on the full tree. Discrete likelihood returns missing-tip state probabilities; parsimony returns optimal states. PIC uses the observed subtree, and PGLS fits complete cases across the response and selected predictors. Neither substitutes inferred means as exact regression observations. UPGMA still requires complete numeric features.

`--traits value,second` produces separate files such as `columns.trait1.nodes.tsv` and `columns.trait2.nodes.tsv`, plus `columns.traits.tsv` mapping trait names to numeric output prefixes and statuses. Each trait has its own observed set and parameters; this is independent batching. CLI batches continue after a column fails and return a nonzero exit code if any fail. Consult the manifest's `status`: prior files for a failed column can still exist. A joint multivariate evolutionary covariance remains future work. `--predictors` supports several PGLS predictors or additional UPGMA features, with singular `--predictor` retained for compatibility.

The [missing-data and uncertainty guide](docs/missing-data-and-uncertainty.md) explains the likelihood, observation flags, API layouts and assumptions. C/Python APIs also accept uncertain discrete observations as likelihood vectors. Continuous measurement-error models, replicate data and uncertain predictors are documented design options for later implementation.

Tree preparation and inference are separate commands, and the two preparation flags also work before an analysis:

```sh
build/obuina prepare --tree examples/tree.nwk --midpoint --prefix results/rooted
build/obuina prepare --tree examples/tree.nwk --extend-tips --prefix results/equal-depth
build/obuina upgma --data examples/traits.tsv --trait value --predictor predictor --prefix results/clustering
```

Scale clustering features yourself when their units differ. A tree built from the response can be useful descriptively, but using that response again for a phylogenetic signal test introduces circularity. Terminal extension holds internal lengths fixed and increases terminal lengths to equal depths; a calibrated dating procedure remains future work.

## Shared outputs

With `--prefix PATH`, analyses produce:

- `PATH.nodes.tsv`: common node IDs, parent IDs, original labels, tip flags, incoming branch lengths and `observed`/`used_for_fit` flags; Brownian means/variances, discrete marginal modes/probabilities or parsimony assignments/state sets where relevant.
- `PATH.summary.tsv`: method, version, trait, missing policy, observed/missing/used tip counts, tree checks, preparation choices and fitted scores/parameters; `state_i` records the class alphabet.
- `PATH.tree.nwk`: the analysed tree, including explicitly requested transformations.
- `PATH.contrasts.tsv` for PIC, or `PATH.coefficients.tsv` for PGLS.

The discrete columns `p_i` and parsimony's comma-separated `optimal_states` use the alphabet in the summary. Independently selected marginal modes may describe an incoherent joint history; the parsimony `state` column is a coherent optimum. Empty fields mean an output does not apply. PIC and PGLS retain the common node table for correspondence and place their numerical results in the additional table. Node IDs refer to the analysed in-memory tree; midpoint rooting and reparsing Newick can change them. The node table's parent IDs and lengths preserve the complete correspondence. Parsimony sets indicate possible optimal states without assigning probabilities. Branchless parsimony reports unavailable depth statistics as `nan` and ultrametricity as `NA`.

Without `--prefix`, the main tables go to standard output and summary metrics to standard error. Prefer a prefix for PIC and PGLS because they produce multiple tables. Existing output files are replaced. Parent output directories must already exist. `prepare` and `upgma` output the tree and summary only.

## C and Python APIs

The public interface is [include/obuina.h](include/obuina.h), with ownership and array conventions documented beside each function. Trait arrays use curupixa node IDs; internal entries are ignored. Result objects own their arrays, errors return `crpx_status`, and constructors clear result pointers on failure. The [C example](examples/basic.c) demonstrates a complete analysis and cleanup.

An installed CMake consumer can use:

```cmake
find_package(obuina CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE obuina::obuina)
```

The [Python binding](python/obuina.py) calls the C ABI directly through standard-library `ctypes`; numerical work runs in C. Results are copied into Python dictionaries before C workspaces are freed. C, Python and the CLI accept several PGLS predictors. `bt_brownian_fit_many` and Python `brownian_many` provide independent multi-column fits; `bt_discrete_asr_likelihoods` and Python `discrete_likelihoods` accept observation-likelihood vectors.

```sh
PYTHONPATH=python python3 examples/python_example.py
PYTHONPATH=python python3 examples/missing_and_uncertain.py
```

```python
from obuina import Tree

with Tree("((A:1,B:2):1,C:2);") as tree:
  fit = tree.brownian({"A": 1.0, "B": 3.0, "C": 4.0})
  print(fit["root_mean"], fit["nodes"])
```

Set `OBUINA_LIBRARY` to the shared library's absolute path when using another build directory or an installation. Manage `Tree` with a context manager. Use a separate tree/analysis per concurrent worker or synchronise mutation and destruction; `ctypes` releases the Python GIL during C calls. Copying owning Python `Tree` objects is blocked; construct `Tree(tree.newick)` for an independent copy.

## Validation and performance

Tests cover exhaustive discrete likelihood/marginal and parsimony enumeration, independent dense Gaussian checks, closed-form GLS, multiple predictors, non-ultrametric and polytomous trees, rooting distance preservation, input errors, CLI outputs and Python ABI parity. See [docs/validation.md](docs/validation.md) for the completed checks and reproducible commands.

Brownian ASR and PIC take linear time/storage in the number of nodes, including observed-subtree handling. The discrete reconstruction passes take $O(NK)$ time/storage for $N$ nodes and $K$ states; bounded rate fitting reuses a prepared upstream binding. UPGMA needs $O(T^2F)$ distance work and $O(T^2)$ storage for $T$ taxa and $F$ features. PGLS currently uses a dense $O(T^3)$ solver on retained complete cases. Midpoint rebuilding uses the upstream tail-aware builder; named-tip validation adds sorting costs recorded in [limitations](docs/limitations.md).

## References and licence

The design draws on the [PastML methods](https://pastml.pasteur.fr/help), [castor reference manual](https://cran.r-project.org/web/packages/castor/castor.pdf), [phangorn ancestral reconstruction interface](https://klausvigo.github.io/phangorn/reference/ancestral.pml.html) and [phytools](https://liamrevell.r-universe.dev/phytools/doc/manual.html). The suggested [prepR4pcm paper](https://peercommunityjournal.org/articles/10.24072/pcjournal.809/) concerns preparation of trees and comparative data; related [package guidance](https://itchyshin.github.io/prepR4pcm/reference/pr_get_tree.html) helped frame explicit tree checks. [Methods and references](docs/methods.md) distinguish the implemented models from these broader packages.

SPDX licence: **GPL-3.0-or-later**, consistent with the linked curupixa dependency. See [LICENSE](LICENSE).
