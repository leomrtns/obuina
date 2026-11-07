---
tags: [phylogenetics, ASR, PCM, references]
---
# Tutorial sources and implementation map

[Tutorial contents](asr-pcm-tutorial.md)

The tutorial's implementation-specific claims were checked against the source and examples below. Equations and worked values explain those calculations; references supply historical context, alternative conventions and extensions. References were checked on 7 November 2026. Links to other software do not imply numerical equivalence to its defaults.

## Obuina implementation

| Topic | Implementation or contract | Existing validation |
|---|---|---|
| Public functions, ownership and array conventions | [obuina.h](../include/obuina.h) | [Interface tests](../tests/test_interfaces.py) |
| Brownian/PIC upstream delegation, missing subtrees and dense PGLS | [continuous.c](../src/continuous.c), upstream [continuous_traits.h](../../curupixa/lib/continuous_traits.h) | [Continuous tests](../tests/test_continuous.c) |
| ER rate fitting and reconstruction wrappers | [discrete.c](../src/discrete.c), upstream [reconstruction.h](../../curupixa/lib/reconstruction.h) and [likelihood.h](../../curupixa/lib/likelihood.h) | [Discrete enumeration tests](../tests/test_discrete.c) |
| Rooting, terminal extension and UPGMA | [trees.c](../src/trees.c) | [Tree tests](../tests/test_trees.c) |
| CLI columns and options | [cli.c](../src/cli.c) | [Interface tests](../tests/test_interfaces.py) |
| Python C-ABI binding | [obuina.py](../python/obuina.py) | [Python example](../examples/python_example.py) |

See the [validation record](validation.md) for the completed numerical and interface checks, and [curupixa findings](curupixa-issues.md) for the dependency version and extension opportunities. The tutorial's examples are small teaching calculations and cannot establish adequacy for a particular empirical application.

The initial tutorial's 11 analysis/preparation CLI commands and three Python examples ran successfully against 0.1.0, with reported values checked against the rounded tables and hand-worked examples. Those complete-data examples remain useful for numerical parity after upstream delegation. New missing-data, independent-batch and uncertain-discrete examples are in the [0.2.0 guide](missing-data-and-uncertainty.md) and [Python example](../examples/missing_and_uncertain.py); the [validation record](validation.md) separates the earlier checks from the adoption checks. Local links and Obsidian maths delimiters can be checked without running the numerical models.

## Foundational work and model references

- [Felsenstein (1985), *Phylogenies and the comparative method*](https://doi.org/10.1086/284325): the foundational independent-contrasts paper, also cited in ape's PIC documentation.
- [Slater (2014), correction and note on comparative palaeontological data](https://doi.org/10.1111/2041-210X.12201): the source of the ultrametric-shortcut warning discussed in the extensions chapter.
- [Bastide et al. (2023), *A Phylogenetic Framework to Simulate Synthetic Interspecies RNA-Seq Data*](https://pmc.ncbi.nlm.nih.gov/articles/PMC11249980/): an explicit fixed-root OU covariance within a broader trait/observation framework.
- [Mitov et al., *Fast likelihood calculation for multivariate Gaussian phylogenetic models with shifts*](https://venelin.github.io/publications/papers/mitov-tpb/): Gaussian transition-based likelihood calculation, including unequal-depth and multifurcating trees.

## Package documentation used for conventions

- [ape: Brownian covariance](https://search.r-project.org/CRAN/refmans/ape/html/corBrownian.html), [Pagel's $\lambda$](https://search.r-project.org/CRAN/refmans/ape/html/corPagel.html), [PIC](https://search.r-project.org/CRAN/refmans/ape/html/pic.html), and [ancestral estimation](https://search.r-project.org/CRAN/refmans/ape/html/ace.html).
- [PastML methods guide](https://pastml.pasteur.fr/help): marginal, joint and ambiguity-retaining discrete reconstruction summaries.
- [phangorn ancestral reconstruction](https://klausvigo.github.io/phangorn/reference/ancestral.pml.html): parsimony ambiguity and ancestral output conventions.
- [phytools Mk model documentation](https://github.com/liamrevell/phytools/blob/master/man/fitMk.Rd): discrete model families and rate parameterisations.
- [PCMBase introduction](https://venelin.github.io/PCMBase/articles/PCMBase.html): Gaussian evolutionary transitions and root/observation specification.
- [geiger reference manual](https://cran.r-project.org/web/packages/geiger/refman/geiger.html): continuous model options and directional-trend identifiability.

The Obuina chapters provide worked derivations from their stated models. They do not assume that options called `marginal`, `joint`, `rate` or `variance` have identical meanings across these packages. Before an external numerical comparison, match the state alphabet, rate normalisation, root assumptions, ML/REML convention, covariance versus correlation and the treatment of polytomies and observations.
