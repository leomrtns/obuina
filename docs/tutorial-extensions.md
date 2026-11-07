---
tags: [phylogenetics, PCM, OU, model-assumptions, tutorial]
---
# Extensions, general trees and model-specific pitfalls

[Tutorial contents](asr-pcm-tutorial.md) · [Previous: discrete traits](tutorial-discrete.md)

**The models in this chapter are theoretical extensions; Obuina 0.2.0 has no OU, Pagel's $\lambda$, directional-trend, early-burst or ordered-CTMC fitting command.** They explain how the implemented methods fit into a broader PCM framework and which details matter when adding models. Their implementation tasks are recorded in [limitations.md](limitations.md). Missing-tip inference and discrete observation-likelihood vectors are now implemented separately; see the [missing-data and uncertainty guide](missing-data-and-uncertainty.md).

## 1. Ornstein–Uhlenbeck evolution

### Branch transition and interpretation

An OU process adds attraction towards a value $\theta$ to continuous random variation:

$$
dY_t=\alpha(\theta-Y_t)\,dt+\sigma\,dW_t,
\qquad\alpha>0.
$$

$\theta$ is the attracting level, $\alpha$ controls how quickly expected states approach it, and $\sigma^2$ controls diffusion. Over an edge of length $t$,

$$
Y_c\mid Y_p=x\sim\mathcal N\left(
\theta+e^{-\alpha t}(x-\theta),\;
\frac{\sigma^2}{2\alpha}(1-e^{-2\alpha t})\right).
$$

The phylogenetic half-life $\log2/\alpha$ is the time required to halve a deviation in the conditional mean. In a non-genomic application, attraction might represent stabilisation around an operating level; interpreting it as biological adaptation requires additional justification. Fitting attraction alone does not identify its causal mechanism.

OU remains Gaussian, so it can use pruning rather than a dense tip covariance. Its branch message has the form $Y_c=a_eY_p+b_e+\epsilon_e$, whereas Brownian motion has $a_e=1,b_e=0$. Extending the Gaussian transitions and root handling can support this calculation on unequal-depth trees. The [PCMBase framework](https://venelin.github.io/PCMBase/articles/PCMBase.html) implements Gaussian transitions of this broader form, and the [authors' algorithm description](https://venelin.github.io/publications/papers/mitov-tpb/) explicitly includes non-ultrametric trees and polytomies.

### Root treatment changes the covariance

Take a fixed root value $x_0$. With tip depths $t_i,t_j$ and shared depth $s_{ij}$,

$$
E(Y_i\mid x_0)=\theta+e^{-\alpha t_i}(x_0-\theta),
$$

$$
\operatorname{Cov}(Y_i,Y_j\mid x_0)
=\frac{\sigma^2}{2\alpha}
e^{-\alpha(t_i+t_j-2s_{ij})}
\left(1-e^{-2\alpha s_{ij}}\right).
$$

One way to derive the covariance is to accumulate shared noise until the MRCA, then attenuate that inherited variation independently along both descendant paths. Equivalently,

$$
\operatorname{Cov}(Y_i,Y_j\mid x_0)
=\sigma^2\int_0^{s_{ij}}
e^{-\alpha(t_i-u)}e^{-\alpha(t_j-u)}\,du.
$$

The fixed-root covariance is also given in the model description of [Bastide et al.'s interspecies simulation framework](https://pmc.ncbi.nlm.nih.gov/articles/PMC11249980/). The integral here explains why both the shared depth and the later path lengths matter.

If instead the root is random with variance $v_0$, add the inherited-root contribution

$$
v_0e^{-\alpha(t_i+t_j)}
$$

to the fixed-root covariance, assuming root variation is independent of subsequent noise. A stationary root has mean $\theta$ and variance $v_0=\sigma^2/(2\alpha)$, giving

$$
\operatorname{Cov}_{\mathrm{stationary}}(Y_i,Y_j)
=\frac{\sigma^2}{2\alpha}e^{-\alpha d_{ij}},
\qquad d_{ij}=t_i+t_j-2s_{ij}.
$$

Stationary-root and fixed-root models are different assumptions about the process before the observed tree begins. An implementation must specify which it fits. With a stationary root, every tip has the same marginal variance despite unequal depths. With a fixed root, tip variance is $\sigma^2(1-e^{-2\alpha t_i})/(2\alpha)$ and generally differs across tips.

As $\alpha\to0$ with fixed root and finite $\sigma^2$, the branch transition and fixed-root covariance approach Brownian motion. The stationary root variance diverges in that limit, so a stationary-root OU likelihood does not approach the fixed-root Brownian likelihood by simply substituting $\alpha=0$. Weak attraction, root uncertainty and limited tree depth can also make parameter estimation difficult.

### Slater's correction: an ultrametric shortcut used on unequal depths

[Slater (2014)](https://doi.org/10.1111/2041-210X.12201) corrected equations and an implementation for models combining OU and Brownian evolution on palaeontological trees. The issue concerned branch-rescaling shortcuts that assumed contemporaneous tips. On an ultrametric tree, $t_i=t_j=T$, so the descendant path distance becomes $2(T-s_{ij})$. Substituting that expression for the actual distance on a non-ultrametric tree gives the wrong OU covariance. The correction concerns that calculation, not a requirement that the OU process itself have equal tip depths.

For a general tree, retain the actual $t_i+t_j-2s_{ij}$. Implement OU using its edge transitions or the appropriate covariance and root model. The original article also addresses mixed-process corrections; those need their own derivations when OU and Brownian regimes change along paths.

### A numerical illustration of the failure

Suppose A and B form a descendant clade whose MRCA with C lies at shared depth 1. Let $t_A=2$, $t_B=3$, $t_C=4$, $\alpha=0.5$ and $\sigma^2=1$. The fixed-root formula above gives

$$
V_{AC}=e^{-2}(1-e^{-1})\approx0.085548,
\qquad
V_{BC}=e^{-2.5}(1-e^{-1})\approx0.051888.
$$

Although $s_{AC}=s_{BC}$, the inherited variation is attenuated for longer along B's path. An ordinary Brownian covariance calculated solely from rescaled edge lengths on that unchanged topology would force those two covariances to be equal, because both pairs have the same shared path. Such rescaling alone cannot reproduce this OU covariance. General Gaussian pruning or direct covariance calculation can.

This illustration uses the transition/covariance equations above and chosen numerical values; it is not a reproduction of Slater's empirical example. It provides a small acceptance case for a future Obuina OU implementation.

## 2. Pagel's $\lambda$ and the covariance/correlation distinction

A useful extension of Brownian residual covariance is

$$
C_\lambda=\lambda C+(1-\lambda)\operatorname{diag}(C),
\qquad0\leq\lambda\leq1.
$$

It retains each tip's variance factor and multiplies each off-diagonal covariance by $\lambda$, matching [ape's definition](https://search.r-project.org/CRAN/refmans/ape/html/corPagel.html). The interval above is a conservative range yielding a positive-definite convex combination when $C$ is positive definite. Some implementations consider larger, tree-dependent admissible ranges, which require separate checking.

At $\lambda=1$, this is the original Brownian covariance. At $\lambda=0$, residuals are independent but retain variances proportional to $t_i$. On a non-ultrametric tree these variances differ, so $\lambda=0$ gives weighted least squares under this formulation. Ordinary least squares follows only if those variances are equal or a different variance model is specified.

This leads to a common cross-package comparison issue. Define $D=\operatorname{diag}(\sqrt{C_{11}},\ldots,\sqrt{C_{nn}})$ and the correlation matrix

$$
R=D^{-1}CD^{-1},\qquad C=DRD.
$$

Using $R$ with a common residual variance fits $\sigma^2R$. Obuina' Brownian PGLS uses $\sigma^2C$. They are equivalent up to a common scale on an ultrametric tree, but differ when diagonal depths vary. A software interface that accepts only a phylogenetic correlation structure needs the corresponding variance weights to reproduce the full Brownian covariance. Check whether an API returns covariance or correlation before comparing coefficients and standard errors; [ape's Brownian interface](https://search.r-project.org/CRAN/refmans/ape/html/corBrownian.html) exposes both.

Estimating $\lambda$ would let the data adjust residual phylogenetic dependence within this covariance family. It would not resolve an incorrect conditional-mean model, circularly inferred tree or inappropriate observation distribution.

## 3. Directional trends, changing rates and regimes

A Brownian process with constant directional drift $b$ has

$$
Y_c\mid Y_p=x\sim\mathcal N(x+bt_e,\sigma^2t_e),
\qquad E(Y_i)=\mu+bt_i.
$$

When all tips have depth $T$, the mean contains only $\mu+bT$. With an unknown root mean and no additional ancestral information, those two contributions cannot be separated from tip means. Unequal meaningful observation times can provide information about the drift. The [geiger manual](https://cran.r-project.org/web/packages/geiger/refman/geiger.html) notes this identifiability restriction for its trend model. Unequal substitution lengths alone do not demonstrate sampling-time differences; the drift's units and interpretation must follow the chosen axis.

A changing diffusion rate can instead preserve a constant mean. For an early-burst/accelerating model with rate $\sigma^2(u)=\sigma_0^2e^{au}$ at root-based time $u$, an edge spanning $u_0$ to $u_1$ has innovation variance

$$
q_e=\int_{u_0}^{u_1}\sigma_0^2e^{au}\,du
=\begin{cases}
\sigma_0^2(e^{au_1}-e^{au_0})/a,&a\ne0,\\
\sigma_0^2(u_1-u_0),&a=0.
\end{cases}
$$

Gaussian pruning can use these edge variance contributions on unequal-depth trees. Different rate regimes can similarly specify edge-dependent innovation variances. These models change accumulation of noise; OU additionally attenuates inherited states through $a_e=e^{-\alpha t_e}$. Their parameters and mechanistic interpretations therefore differ. Fitting regime boundaries, rate parameters or transitions between processes introduces additional inference beyond the currently fixed Brownian model.

## 4. Ordered discrete likelihood

An adjacent-state CTMC gives finite integer levels a probabilistic ordering. One possible generator is

$$
Q_{i,i+1}=u_i,\qquad Q_{i,i-1}=d_i,\qquad
Q_{ij}=0\ \text{for }|i-j|>1,
$$

with diagonals chosen so each row sums to zero and only existing boundary neighbours included. Its transition matrix is still $P(t)=e^{Qt}$. Even when $Q_{0,2}=0$, a positive-length edge can have $P_{0,2}(t)>0$ through successive changes.

The waiting times and directional rates specify more than ordered parsimony's step costs. Uniform stationary frequencies are a particular outcome of suitable symmetric rates, rather than an automatic property of every ordered generator. A bounded alphabet also imposes endpoint behaviour; an unbounded count process requires a different formulation or a carefully checked truncation. [phytools' ordered model documentation](https://github.com/liamrevell/phytools/blob/master/man/fitMk.Rd) illustrates additional rate-parameterisation choices.

General generators also need a reliable transition calculation. A dense matrix exponential can increase cost relative to the ER shortcut. A future implementation should verify stochastic row sums, nonnegative transition probabilities, zero-time identity, long-time limits, root treatment and tiny-tree likelihood/marginal results against independent calculations.

## 5. Non-Gaussian regression and observations

For binary, count or ordinal responses, a phylogenetic latent-variable model is often more natural than Gaussian residuals on the observed scale. A possible structure is

$$
\mathbf u\sim\mathcal N(0,\tau^2C),\qquad
g\{E(Y_i\mid\mathbf u)\}=\mathbf x_i^T\boldsymbol\beta+u_i,
$$

with a specified observation distribution and link $g$. Examples include a Bernoulli response with a logit link or a Poisson response with a log link. Such a model conditions observations on latent phylogenetic effects. It generally requires integration or approximation over those effects, and its observed-scale covariance is not simply $\tau^2C$.

This addresses regression of a response on predictors. A discrete CTMC instead describes state changes along the tree. Although both can analyse binary tip observations, they encode different hypotheses and yield different ancestral/regression interpretations. Neither is fitted by the current `pgls` command. Extra-Poisson variation, excess zeros, observation error and ordinal threshold identifiability would require further choices.

For replicated measurements, distinguish evolution of a latent taxon mean from sampling variation around it. Treating each noisy continuous measurement as an exact tip state can distort evolutionary rates and certainty. Obuina' current exact-observation assumption for supplied continuous values should be reported alongside results; unknown continuous values are handled separately. Discrete observation-likelihood vectors can already represent an explicit noisy classification model. The [uncertainty guide](missing-data-and-uncertainty.md#8-continuous-measurement-uncertainty-proposed-next-step) develops continuous error and future API options.

## 6. Apply the models to non-genomic data

A tree-based model assumes that shared paths represent a source of shared variation. That can be useful beyond genetic sequences when branching inheritance is plausible. Horizontal transfer, recombination, cultural exchange or convergent environments can create dependence the single tree does not describe. Independent predictors or observation models may account for some of that structure; network or other covariance models may be needed for the rest.

The current UPGMA command can construct a descriptive dendrogram from Euclidean metadata distances:

```sh
build/obuina upgma --data examples/traits.tsv --trait value --predictor predictor --prefix results/tutorial/clustering
```

Use this only after the [tutorial output-directory setup](asr-pcm-tutorial.md#4-set-up-the-worked-examples). Features with different units need a justified scaling. If the response helps construct its own tree, that tree cannot serve as independent evidence for the response's phylogenetic signal or provide ordinary fixed-covariance regression calibration without accounting for that construction. An independent tree or explicitly joint inference would address a different statistical design.

## 7. A reporting template

For a reproducible analysis, record the following genuinely separate choices:

| Choice | Information to retain |
|---|---|
| Dependence structure | Tree source, rooting justification, branch units, unequal depths and any transformations |
| Observations | Trait units/coding, state alphabet, missingness, transformations and measurement assumptions |
| Process | Brownian or ER likelihood, ordered/unordered costs, root treatment and fixed/fitted rates |
| Calculation | Software version, CLI/API invocation, estimator conventions, numerical/boundary warnings |
| Output interpretation | Conditional ancestral uncertainty, regression coefficients or parsimony possibilities; sources of uncertainty excluded |

Preserve the input tree and metadata alongside output prefixes. The common node table gives the topology/length correspondence; state mappings in summaries make categorical results interpretable. Future model additions should expand the process declaration and diagnostics while retaining that correspondence.
