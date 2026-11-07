# Methods and assumptions

For step-by-step derivations, comparisons and runnable examples, read the [ASR and PCM tutorial](asr-pcm-tutorial.md) and [missing-data and uncertainty guide](missing-data-and-uncertainty.md). This page records Obuina 0.2.0's method/API conventions.

Obuina treats metadata as observed tip traits. Numerical classes need an explicit evolutionary interpretation: a count can be modelled as continuous under Brownian motion, as an ordered discrete state under parsimony, or eventually through a count process. These choices impose different assumptions even when they use the same observations.

## Trees, rooting and branch lengths

Newick supplies a computational root, including for a nominally unrooted tree. Users must justify its biological interpretation. `--midpoint` roots at the midpoint of the longest weighted tip-to-tip path. It preserves taxon distances and assumes that path lengths provide useful information about root placement. Rate variation, asymmetrical sampling and long branches can undermine that assumption. Outgroup and dated rooting, and minimum-variance alternatives, remain on the [work queue](limitations.md).

Midpoint rooting suppresses unlabelled degree-two vertices and removes unary stems above the first branching vertex, including their labels. It retains labelled internal vertices elsewhere. A retained labelled vertex may become unary, which PIC rejects. Retained internal labels stay attached to vertices; users should reassess labels that represent rooted-clade support after rerooting. An all-zero tree has no unique informative midpoint and is returned as a clone. A root stem length is ignored by analyses. Trees may contain hard polytomies; PIC requires exactly two children at every internal node.

Ultrametricity means equal summed root-to-tip lengths within a relative numerical tolerance. Positive-length Brownian ASR, PIC and Brownian PGLS permit unequal tip depths; the equal-rates CTMC also permits unequal depths. Interpreting lengths as time requires a defensible temporal scale, particularly for serially sampled or fossil tips. Other future PCM models may impose ultrametricity or calibration constraints.

`bt_tree_extend_tips` solves a restricted linear minimisation. With fixed internal paths and original tip depths $d_i$, choose terminal increments $a_i\geq0$ and a common depth $H$:

$$
\min \sum_i a_i,\qquad d_i+a_i=H.
$$

The solution is $H=\max_i d_i$ and $a_i=H-d_i$. It changes distances between tips. General branch-length fitting, calibrated dating and rate smoothing need separate solvers. The [phytools ultrametric guidance](https://liamrevell.r-universe.dev/phytools/doc/manual.html) likewise distinguishes simple coercion from formal dating.

## Brownian model and continuous ASR

On an edge of length $t$, Brownian change is normally distributed with mean zero and variance $\sigma^2t$. Let $C_{ij}$ be the shared root-to-MRCA length, with $C_{ii}$ equal to the root-to-tip depth. The tip model is

$$
\mathbf y\sim\mathcal N(\mu\mathbf1,\sigma^2 C).
$$

For the observed tips, let $T$ be their count, $\mathbf y$ their data and $C$ their covariance-factor submatrix relative to the original tree root. The fitted root mean and diffusion variance are

$$
\hat\mu=\frac{\mathbf1^TC^{-1}\mathbf y}{\mathbf1^TC^{-1}\mathbf1},\qquad
\hat\sigma^2=\frac{(\mathbf y-\hat\mu\mathbf1)^TC^{-1}(\mathbf y-\hat\mu\mathbf1)}{T}.
$$

The reported profile log likelihood includes $\log|C|$ and all Gaussian constants. Curupixa's Gaussian pruning evaluates it in $O(N)$ time without constructing $C$. Obuina maps missing observations (`NAN` in C) to an observed-descendant tree, retaining the original root and relevant lengths, then delegates fitting/reconstruction upstream. It maps results back to original IDs and predicts entirely unobserved subtrees by propagating the parent's conditional mean and adding $\hat\sigma^2t$ per edge to its variance. The fit's denominator counts observed tips only.

A downward and an upward pass compute every retained node's conditional mean and variance using a flat root prior and the fitted $\hat\sigma^2$. Thus internal estimates use all observed tips, including observations outside the node's descendant subtree. Exactly observed tip variances are zero; missing tip predictions have conditional process uncertainty.

Conditional ASR variances include uncertainty in the ancestral root state under its flat prior. They exclude uncertainty from estimating the diffusion variance, the tree, measurement and the evolutionary model. A normal interval `mean ± 1.96*sqrt(variance)` is conditional on those assumptions. Positive edge lengths, at least two observed tips with distinct finite values and positive fitted variance are required. All-missing, one-observation and constant-observation fits return an error. Infinities are invalid, distinct from missing `NAN`. The C result's `tips` counts observed tips and `nodes` counts original nodes. No trend, OU attraction, continuous measurement error or joint multivariate diffusion is fitted.

## Independent contrasts

At a bifurcation, with descendant estimates $x_a,x_b$ and their accumulated variances $v_a,v_b$, PIC reports $(x_a-x_b)/\sqrt{v_a+v_b}$. The propagated estimate uses inverse-variance weights and propagated variance $v_av_b/(v_a+v_b)$; the incoming edge is then added. The output `variance` is $v_a+v_b$ before standardisation. The sign follows child order and can reverse after equivalent Newick reordering. Regression of contrasts conventionally uses the origin because the transformation removes an intercept. Obuina delegates the contrast kernel to curupixa. Missing tips are omitted, unary retained paths are contracted and contrasts map back to original branching-node IDs; no reconstructed means are inserted as observations. The original tree must be binary with positive lengths. The present CLI exports contrasts for downstream analysis.

## Brownian PGLS

PGLS fits $\mathbf y=X\boldsymbol\beta+\boldsymbol\varepsilon$, with $\boldsymbol\varepsilon\sim\mathcal N(0,\sigma^2C)$ and the same rooted covariance definition. An intercept is added automatically. A tip with missing (`NAN`) response or any selected predictor is omitted. $T$ here counts complete cases and $C$ is their submatrix relative to the original root. Cholesky whitening followed by reorthogonalised QR gives the coefficients. The reported variance is ML residual SSE/$T$; coefficient standard errors use SSE/$(T-p)$, where $p$ includes the intercept. No hypothesis-test probabilities or model selection criterion is currently emitted.

The covariance is dense. Cost is $O(T^3)$ time and $O(T^2)$ memory, with additional $O(Tp)$ design storage. Full-rank predictors, $T>p$, finite retained data and positive residual variance are required. Numerically exact zero-residual fits return an error; round-off can leave tiny positive residuals in almost exact fits. No single imputation is performed. The C result's `tips` counts used cases. Pagel's $\lambda$, alternative covariance models, measurement error, joint missing-predictor inference and tree-pruning GLS acceleration remain future work.

## Equal-rates discrete likelihood and ASR

For $K\geq2$ states and multiplier $r$, the generator has $Q_{ij}=r/(K-1)$ for $i\ne j$ and $Q_{ii}=-r$. This uses curupixa's unit expected exit-rate normalisation. For an edge of length $t$:

$$
P_{ij}(t)=\frac1K+\left(\mathbb1_{i=j}-\frac1K\right)\exp\left(-\frac{Krt}{K-1}\right).
$$

The root prior is uniform. Curupixa's reconstruction API computes log-space likelihoods and full all-tip marginal probabilities. Obuina wraps those owned results and supplies the equal-rates fitting policy. This is marginal reconstruction with plug-in parameters. Maximising each node's marginal independently can give a different scenario from maximising the joint probability of all ancestral states. [PastML](https://pastml.pasteur.fr/help) explains the distinction between marginal, joint and MPPA reconstructions; this version returns full marginals and their individual modes.

Unknown observations permit every state. Nonnegative lengths, zero-length transitions and polytomies are accepted. Contradictory observations linked only by zero-length paths have likelihood zero; conditioning then has undefined marginals and returns an error. A completely unknown trait has uniform marginals and does not identify a rate. All-zero edge lengths likewise cannot identify a rate.

Optional rate fitting evaluates a 65-point logarithmic grid and refines the best neighbouring interval with golden-section search. Bounds are $10^{-8}/\bar t$ and $100/\bar t$, where $\bar t$ is the arithmetic mean over all non-root edges, including zero edges. One prepared curupixa likelihood binding avoids repeated mapping/traversal setup during the search. Estimates on a boundary or likelihood plateaux require interpretation; a scalar grid/refinement does not guarantee discovery of every narrow global optimum. The fixed-rate API avoids the repeated evaluations. The reported probabilities exclude fitted-rate uncertainty. Natural-number labels do not impose an ordered CTMC; that generator remains deferred.

### Uncertain discrete observations

`bt_discrete_asr_likelihoods` supplies finite nonnegative weights $w_i(s)=p(d_i\mid Z_i=s)$ in node-major `[node][state]` layout. Internal rows are ignored and each tip requires at least one positive weight. Pruning starts from that vector. One-hot vectors represent exact observations, all ones represent unknown data, and several positive entries can represent ambiguity or a specified noisy observation model. These are observation likelihoods rather than automatically calibrated classifier posteriors. Scaling a whole vector by a positive constant preserves marginals and rate estimates but changes the log-likelihood offset. The Python equivalent is `discrete_likelihoods`; the TSV interface currently exposes exact/unknown states only. See the [uncertainty guide](missing-data-and-uncertainty.md#6-discrete-observation-uncertainty-is-available-through-the-api).

## Parsimony

Sankoff dynamic programming minimises the summed transition costs over all edges. Unordered transitions cost zero when equal and one otherwise. Ordered transitions cost $|i-j|$ using the declared state indices. An observed tip restricts its allowed state; an unknown tip allows every state. The implementation supports hard polytomies, unary nodes and a one-state alphabet.

Curupixa's reconstruction API supplies linear min-plus transforms, making both the unordered and ordered passes $O(NK)$ rather than using a dense $K\times K$ cost matrix. An upward pass marks every state occurring at a node in at least one global optimum. These sets may include choices that cannot all occur together. A separate traceback returns one globally coherent optimum, resolving ties by the lowest allowed index. Obuina maps these results to its public arrays. The state sets indicate possibilities without weighting optimal histories. Branch lengths, including absent lengths, do not enter the cost.

## Missing input policy and multi-column analyses

The CLI defaults to strict metadata coverage. `--missing estimate` permits absent tree-tip rows and continuous `?`, `NA` or empty values; extra or duplicate names remain errors. Explicit discrete unknown markers are accepted with either policy. Python uses `allow_missing=True` for omitted keys and continuous `None`; C uses `NAN`/`SIZE_MAX` or all-ones likelihood vectors as specified above. `.nodes.tsv` records observed/used flags and summaries include the selected trait and observed/missing/used counts.

`--traits a,b` runs independent single-trait analyses with separate prefixes and a `.traits.tsv` manifest. `--predictors x,z` selects several PGLS predictors or extra UPGMA features. The original singular options remain supported. `bt_brownian_fit_many` and Python `brownian_many` batch independent columns; the C batch input is node-major `[node][trait]` and owned batch results are returned only when all traits succeed. This does not estimate cross-trait covariance or borrow information between columns. UPGMA still requires complete numeric features. The [missing-data guide](missing-data-and-uncertainty.md) derives the statistical distinctions.

The [phangorn reconstruction documentation](https://klausvigo.github.io/phangorn/reference/ancestral.pml.html) highlights non-uniqueness of parsimony reconstructions; [castor](https://cran.r-project.org/web/packages/castor/castor.pdf) provides a broader range of reconstruction and transition-cost models. Custom costs, ACCTRAN/DELTRAN and scenario counting are deferred here.

## Reference context

[Revell (2024), phytools 2.0](https://doi.org/10.7717/peerj.16505) describes a much broader PCM ecosystem. The suggested [Nakagawa et al. (2026), prepR4pcm](https://peercommunityjournal.org/articles/10.24072/pcjournal.809/) focuses on preparing comparative trees and data; its [tree-retrieval guidance](https://itchyshin.github.io/prepR4pcm/reference/pr_get_tree.html) makes tree checks and transformations explicit. The publisher article endpoint returned HTTP 403 during development, so detailed package behaviour was checked through the official package guidance and publisher index rather than inferred from the unavailable article text.

The Gaussian and contrast calculations follow the standard Brownian comparative framework of Felsenstein (1985), [Phylogenies and the comparative method](https://doi.org/10.1086/284325). All algorithms here operate on a supplied topology after any explicitly requested preparation. They introduce no topology search or Bayesian sampling.
