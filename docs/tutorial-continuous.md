---
tags: [phylogenetics, ASR, PCM, Brownian-motion, PGLS, tutorial]
---
# Continuous traits: Brownian ASR, contrasts and PGLS

[Tutorial contents](asr-pcm-tutorial.md) · [Next: discrete traits](tutorial-discrete.md)

This chapter uses one evolutionary covariance in three ways: estimate the process and reconstruct ancestral values, form independent contrasts, and fit a regression. Obuina implements all three. Throughout, $n$ is the number of observed tips and $N$ the number of all nodes. Branch lengths need not produce equal tip depths.

## 1. Brownian motion along a tree

Under Brownian motion, a child inherits its parent's state plus an independent Gaussian increment:

$$
Y_v=Y_u+\epsilon_e,\qquad
\epsilon_e\sim\mathcal N(0,\sigma^2t_e),
$$

where $e=(u,v)$ is the edge from parent $u$ to child $v$. Different edges have independent increments. The diffusion variance rate $\sigma^2$ has units of squared trait units per branch unit. Conditional on the ancestor, the expected change is zero; longer edges permit more variation.

An entire tip vector therefore has distribution

$$
\mathbf y\mid\mu,\sigma^2,\mathcal T
\sim\mathcal N(\mu\mathbf1,\sigma^2C),
\qquad C_{ij}=s_{ij},\quad C_{ii}=t_i.
$$

Two tips share increments along the path from the root to their MRCA. Their subsequent independent increments increase their individual variances. This construction gives the shared-depth covariance also documented by [ape's Brownian model](https://search.r-project.org/CRAN/refmans/ape/html/corBrownian.html).

### Assumptions and interpretation

The model assumes a single constant diffusion rate, zero directional drift, exact tip observations and independent increments after branching. It permits any real state. A positive fitted variance does not establish that biological evolution involved neutral genetic drift: Brownian motion is also a statistical model for cumulative variation in non-genomic measurements.

Unequal tip depths give unequal marginal variances. Forcing the tree to equal depths would change that model. The implementation requires positive edges and at least two distinct finite observed values; it supports polytomies and unary nodes. Missing observations can be integrated out for fitting and predicted with conditional variances. Use `--missing estimate` in the CLI or `allow_missing=True` in Python. Continuous measurement error and joint multivariate trait diffusion remain deferred; see [missing data and uncertainty](missing-data-and-uncertainty.md).

## 2. Fit the root mean and diffusion variance

With a fixed but unknown root mean $\mu$, define the following on the observed-tip vector and its covariance submatrix when some tips are missing. Here $n$ counts observed tips:

$$
S(\mu)=(\mathbf y-\mu\mathbf1)^TC^{-1}(\mathbf y-\mu\mathbf1).
$$

The Gaussian log likelihood is

$$
\ell(\mu,\sigma^2)
=-\frac12\left[n\log(2\pi)+n\log\sigma^2+\log|C|+
\frac{S(\mu)}{\sigma^2}\right].
$$

Differentiating with respect to $\mu$ and $\sigma^2$ gives

$$
\widehat\mu=\frac{\mathbf1^TC^{-1}\mathbf y}{\mathbf1^TC^{-1}\mathbf1},
\qquad
\widehat\sigma^2_{\mathrm{ML}}=\frac{S(\widehat\mu)}{n}.
$$

The mean weighs observations according to covariance, so a dense group of closely related tips does not count as the same number of independent observations. The rate describes remaining evolutionary variation around this fitted mean.

Obuina reports this **maximum-likelihood** variance, with divisor $n$. Restricted maximum likelihood, REML, accounts for the fitted mean and gives divisor $n-1$ in this intercept-only Brownian setting. Thus package comparisons must match the estimator. The [ape ASR documentation](https://search.r-project.org/CRAN/refmans/ape/html/ace.html) explicitly distinguishes ML and REML; Obuina currently fits ML.

### How Gaussian pruning avoids a dense inverse

At an observed tip, the exact observation is a message with mean $m_i=y_i$ and variance factor $v_i=0$. Moving a child's message across its incoming edge adds that edge length: $V_i=v_i+t_i^{\mathrm{edge}}$. Missing-only subtrees carry no observation information; Obuina delegates fitting to curupixa on the observed-descendant tree, retaining the original root and relevant path lengths. This yields the same observed-tip marginal likelihood without inserting fitted values as observations.

For two child messages at the same parent, multiplication of their Gaussian kernels yields

$$
m=\frac{m_a/V_a+m_b/V_b}{1/V_a+1/V_b},\qquad
v=\frac{V_aV_b}{V_a+V_b}.
$$

Their disagreement contributes

$$
\frac{(m_a-m_b)^2}{V_a+V_b}
$$

to the profiled quadratic $S$. Each merge also contributes $\log(V_a+V_b)$ to the covariance log determinant. At the root, adding $\log v_{\mathrm{root}}$ completes $\log|C|$. Repeating merges accommodates polytomies. Each edge participates in a bounded number of operations, giving $O(N)$ time and storage.

The downward messages contain descendant information. At the root all tips are descendants, so its downward mean is $\widehat\mu$. Other internal downward means generally need an outside-subtree contribution before they become full ASR estimates.

## 3. Reconstruct every ancestral state

For reconstruction, Obuina places a flat prior on the ancestral root state and plugs in $\widehat\sigma^2_{\mathrm{ML}}$. The conditional distribution of all unobserved states $\mathbf z$ (internal nodes and any missing tips) is proportional to

$$
p(\mathbf z\mid\mathbf y,\widehat\sigma^2,\mathcal T)
\propto\exp\left[-\frac{1}{2\widehat\sigma^2}
\sum_{(u,v)\in E}\frac{(x_v-x_u)^2}{t_{uv}}\right],
$$

where $\mathbf x$ contains the unknown states and fixed observed tips. Positive edges and observations make this distribution proper, despite the flat root prior.

Expanding the sum gives an unobserved-node precision matrix $K$. Each incident edge contributes $1/t_e$ to an unobserved node's diagonal; an edge joining two unobserved nodes contributes $-1/t_e$ to the corresponding off-diagonal. For an unobserved node adjacent to an observed tip, that edge contributes $y_i/t_e$ to an information vector $\mathbf h$. Therefore

$$
E(\mathbf z\mid\mathbf y)=K^{-1}\mathbf h,
\qquad
\operatorname{Cov}(\mathbf z\mid\mathbf y)=\widehat\sigma^2K^{-1}.
$$

These equations explain the result. The upstream implementation combines Gaussian messages in downward and upward passes. For a child, the upward message includes its parent's outside information and all its siblings; it excludes that child's own downward contribution. Multiplying the two messages gives its full conditional mean and variance. Prefix/suffix products keep this calculation linear even at polytomies. Obuina maps retained results to original IDs, then propagates means and increasing variances down missing-only subtrees.

### Fitting and reconstruction use different root operations

The reported `log_likelihood` profiles the tip likelihood over a fixed root parameter $\mu$, as in section 2. The ASR distributions then integrate the ancestral root state under a flat prior. Consequently the reported root mean agrees with $\widehat\mu$, while the root has positive conditional variance. Fixing the root at $\widehat\mu$ during reconstruction would give different uncertainty.

The ASR variance includes uncertainty about the root state and other ancestral states. It excludes uncertainty in estimating $\sigma^2$, the tree, branch lengths and measurement. A conditional Gaussian interval is

$$
m_v\pm1.96\sqrt{\operatorname{Var}(Y_v\mid\mathbf y)}.
$$

That interval conditions on the fitted rate. Exactly observed tip variances are zero; missing tip predictions have positive uncertainty. With constant observed data, the ML variance lies on the boundary at zero; this implementation returns an error rather than a finite positive-variance fit. At least two observed tips are required to fit that variance. The [missing-tip derivation](missing-data-and-uncertainty.md#3-brownian-prediction-includes-root-state-uncertainty) gives the corresponding Gaussian conditional formula, including root-state uncertainty.

## 4. A small calculation by hand

Consider

```text
((A:1,B:2)I:1,C:4)R;
```

and observations $y_A=1$, $y_B=4$, $y_C=7$. In tip order A, B, C,

$$
C=\begin{pmatrix}2&1&0\\1&3&0\\0&0&4\end{pmatrix},\qquad
C^{-1}=\begin{pmatrix}3/5&-1/5&0\\-1/5&2/5&0\\0&0&1/4\end{pmatrix}.
$$

The tree is non-ultrametric: depths are $2,3,4$. Calculation gives

$$
\widehat\mu=\frac{59}{17}\approx3.470588,\qquad
S=\frac{126}{17},\qquad
\widehat\sigma^2=\frac{42}{17}\approx2.470588.
$$

Since $|C|=20$, the profiled log likelihood is approximately $-7.111366$. The full ASR results are

| Node | Conditional mean | Conditional variance |
|---|---:|---:|
| R | $59/17\approx3.470588$ | $840/289\approx2.906574$ |
| I | $44/17\approx2.588235$ | $420/289\approx1.453287$ |

Using only A and B gives I a downward mean of $2$ and variance factor $2/3$. Including C through the upward pass shifts its full mean to $44/17$. This difference is why a subtree-only reconstruction differs from an estimate conditioned on all tips.

To reproduce from the project root, run the following Python block with `PYTHONPATH=python` set in the environment:

```python
from obuina import Tree

with Tree("((A:1,B:2)I:1,C:4)R;") as tree:
  fit = tree.brownian({"A": 1.0, "B": 4.0, "C": 7.0})
  print(fit["root_mean"], fit["sigma2"], fit["log_likelihood"])
  for node in fit["nodes"]:
    if not node["is_tip"]:
      print(node["label"], node["mean"], node["variance"])
```

## 5. Independent contrasts

A contrast summarises a difference between two independently diverging descendant lineages. For the same descendant summaries and accumulated variance factors as above,

$$
c_v=\frac{m_a-m_b}{\sqrt{V_a+V_b}}.
$$

The parent summary uses the same inverse-variance mean and variance-factor recurrence as Gaussian pruning. Under the specified Brownian process, the $n-1$ contrasts of a bifurcating tree are independent Gaussian variables with variance $\sigma^2$. The standardisation removes tree-dependent variance factors; it does not divide by the estimated diffusion standard deviation.

For the hand example, the contrast at I is $-3/\sqrt3\approx-1.732051$. At R it is $-5/\sqrt{17/3}\approx-2.100420$. Their squared sum is $126/17$, the same $S$ used in the Brownian fit. The output `variance` is the denominator's squared factor, $3$ or $17/3$, in branch units.

```sh
build/obuina pic --tree examples/tree.nwk --data examples/traits.tsv --prefix results/tutorial/pic
```

The sign follows child order. For a two-trait analysis, calculate both traits' contrasts on the same tree with the same child ordering; swapping a node's children then reverses both signs consistently. Regress contrasts through the origin. The [ape PIC examples](https://search.r-project.org/CRAN/refmans/ape/html/pic.html) illustrate this convention.

PIC's propagated internal summaries use descendants only. Obuina exports contrasts through curupixa's `pic` kernel; use `brownian` for full all-tip ASR. The original input must be binary with positive edges. With missing data, Obuina retains observed tips, contracts paths with one retained descendant and calculates $n_O-1$ contrasts without imputed observations. Original branching-node IDs are preserved in the returned correspondence. Unequal tip depths are supported.

## 6. PGLS: use Brownian covariance for regression residuals

For a response and predictors, specify

$$
\mathbf y=X\boldsymbol\beta+\boldsymbol\varepsilon,\qquad
\boldsymbol\varepsilon\sim\mathcal N(0,\sigma^2C).
$$

$X$ includes an intercept column. The covariance applies to **residuals after the predictors**, so the response itself need not have a constant Brownian mean. Standard regression interpretation also requires a suitable conditional-mean specification, including that unexplained residuals have zero expectation given the predictors.

Generalised least squares minimises

$$
S(\boldsymbol\beta)=(\mathbf y-X\boldsymbol\beta)^TC^{-1}(\mathbf y-X\boldsymbol\beta),
$$

with solution

$$
\widehat{\boldsymbol\beta}=(X^TC^{-1}X)^{-1}X^TC^{-1}\mathbf y.
$$

The formula expresses the estimator. Obuina constructs $C$, calculates a Cholesky factor $C=LL^T$, whitens $\mathbf y^*=L^{-1}\mathbf y$ and $X^*=L^{-1}X$, then solves the whitened regression using reorthogonalised QR. This avoids explicitly forming $C^{-1}$ or solving squared-condition-number normal equations.

For $p$ fitted coefficients including the intercept, Obuina reports

$$
\widehat\sigma^2_{\mathrm{ML}}=S(\widehat{\boldsymbol\beta})/n,
\qquad
\widehat{\operatorname{Cov}}(\widehat{\boldsymbol\beta})=
\frac{S(\widehat{\boldsymbol\beta})}{n-p}(X^TC^{-1}X)^{-1}.
$$

The standard errors use the latter matrix's diagonal. The reported Gaussian log likelihood uses the ML variance, including the determinant and normalising constants. Matching these two variance conventions matters when reproducing output.

The model requires $n>p$, a full-rank retained design, finite retained data, positive lengths and positive fitted residual variance. A row with missing response or any selected predictor is omitted; here $n$ and $C$ refer to complete cases and their original-root covariance submatrix. No inferred mean is treated as an observed regression datum. Numerically exact fits can produce a zero-variance boundary error; nearly exact fits can leave very small positive residuals through round-off. Extreme scales or poor conditioning can also cause a recoverable numerical error. The current dense calculation costs $O(n^3)$ time and $O(n^2)$ covariance storage, unlike the linear Brownian ASR calculation.

### Connection to contrasts and ASR

For the same Brownian covariance on a bifurcating tree, a contrast matrix $A$ satisfies

$$
A\mathbf1=0,\qquad ACA^T=I_{n-1}.
$$

The first property removes the intercept; the second removes phylogenetic covariance. Regressing response contrasts on predictor contrasts through the origin recovers the PGLS slopes for that matching model. The intercept still needs a separate calculation. This equivalence depends on matching lengths, covariance and observation assumptions.

With no predictors, PGLS reduces to fitting the Brownian root mean: its intercept, ML variance and likelihood match `brownian`. Its coefficient standard error uses $S/(n-1)$, while Brownian ASR uses the ML plug-in rate $S/n$ for ancestral variances. With predictors, a PGLS intercept describes the response when the predictors equal zero; it does not automatically estimate an ancestral response. Ancestral predictor values and a model connecting them to the response would be needed for that interpretation.

## 7. Run and interpret the supplied example

After the [tutorial setup](asr-pcm-tutorial.md#4-set-up-the-worked-examples), run

```sh
build/obuina brownian --tree examples/tree.nwk --data examples/traits.tsv --prefix results/tutorial/brownian
build/obuina pgls --tree examples/tree.nwk --data examples/traits.tsv --predictor predictor --prefix results/tutorial/pgls
```

Expected values, rounded from the current implementation:

| Quantity | Value |
|---|---:|
| Brownian fitted root mean | 2.295301 |
| Brownian diffusion variance rate, ML | 0.972354 |
| Brownian tip log likelihood | -7.167472 |
| PGLS intercept | 0.199318 |
| PGLS predictor coefficient | 2.005860 |
| Intercept standard error | 0.316748 |
| Predictor coefficient standard error | 0.226797 |
| PGLS residual diffusion variance, ML | 0.024242 |

The predictor coefficient describes the conditional change in response per predictor unit under the chosen residual covariance. The four-tip example is designed to demonstrate calculations; it provides little information for substantive inference. A small standard error in this example reflects its small fitted residual variation and the assumed model. It does not establish causation or rule out measurement error, a wrong tree or omitted predictors.

The CLI accepts one predictor with `--predictor` or several with `--predictors x,z`, matching the multi-predictor Python/C APIs. For the Python equivalent of the single-predictor example, with `PYTHONPATH=python`:

```python
from obuina import Tree

y = {"A": 1.2, "B": 2.8, "C": 2.1, "D": 4.5}
x = {"A": 0.5, "B": 1.5, "C": 0.9, "D": 2.0}
with Tree("((A:1,B:2)AB:1,(C:0.8,D:1.4)CD:1.5)root;") as tree:
  print(tree.pgls(y, [x]))
```

## 8. Diagnose assumptions before adding complexity

Inspect whether large contrasts concentrate in a clade, vary systematically with their expected standard deviations, or coincide with unusual sampling/measurement. Such patterns may indicate heterogeneous rates, an unsuitable branch scale or observation errors. For PGLS, inspect whitened residuals, fitted relationships and influential taxa; phylogenetic correction does not repair a misspecified mean model.

Compare plausible trait transformations and independent defensible trees where available. Keep the response, its units and its observation model consistent when comparing likelihoods: raw and log-transformed responses need a Jacobian adjustment to compare densities on the same original scale. Do not compare a continuous log likelihood directly with a categorical likelihood or parsimony score. Formal diagnostic, model-comparison and uncertainty-propagation tools remain on the [work queue](limitations.md).
