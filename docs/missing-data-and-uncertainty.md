---
tags: [phylogenetics, ASR, PCM, missing-data, observation-uncertainty, obuina]
---
# Missing data, several metadata columns and uncertainty

Obuina **0.2.0** can estimate unobserved continuous or discrete tip states under the same evolutionary models used for ASR. It can also analyse several columns independently. This page distinguishes missing observations, uncertain observations and uncertainty in fitted evolutionary parameters, and specifies the implemented behaviour separately from proposed extensions.

[Tutorial contents](asr-pcm-tutorial.md) · [Method contracts](methods.md) · [Current limitations](limitations.md)

## 1. What each method does when a tip has no observation

| Method | Implemented treatment | Information used |
|---|---|---|
| Brownian ASR | Fit the observed-tip marginal likelihood; predict every missing tip with a mean and variance | Other observed tips in the same trait column and the supplied tree |
| Discrete ER ASR | Allow every state at the missing tip; return its conditional probability vector | Other observations or observation-likelihood vectors in the same trait |
| Parsimony | Allow every state; return a coherent optimal assignment and global optimal-state sets | Declared costs and observed states |
| PIC | Calculate contrasts on the observed-tip subtree, contracting paths with only one retained descendant | Observed continuous values; no inserted predictions |
| PGLS | Omit tips missing the response or any selected predictor; fit the complete-case covariance submatrix | Retained response/predictor rows and their original-root Brownian covariance |
| UPGMA | Require complete selected feature rows | Euclidean feature distances; no missing-feature distance model |

The model remains defined on the supplied tree. Missing values do not contribute fictitious observations. Replacing a missing value by one fitted mean and treating that mean as exact in a subsequent regression generally understates uncertainty and changes the effective amount of information.

### CLI and Python policy

The CLI keeps strict checking by default. An omitted taxon row is an error unless `--missing estimate` is supplied. Explicit discrete unknown markers (`?`, `NA` or an empty field) remain valid without that flag. Continuous unknown markers require `--missing estimate`. Extra taxa, duplicate rows and duplicate tip names remain errors; the missing-data option does not relax name matching for supplied rows.

In Python, `allow_missing=True` permits omitted taxon keys. For continuous traits it also permits `None` observations. Explicit `None` was already supported as an unknown discrete state; an omitted discrete key needs the flag. C callers use `NAN` for missing continuous values, `SIZE_MAX` for an unknown exact discrete state, or an all-ones discrete observation-likelihood vector. Positive/negative infinity remains invalid.

Missingness is evaluated separately for each selected column. A taxon can inform one trait and remain unobserved for another. In PGLS, a row contributes only when its response and every selected predictor are observed.

## 2. Brownian fitting uses only observed observations

Partition the tips into observed indices $O$ and missing indices $M$. Under the Brownian model,

$$
\begin{pmatrix}\mathbf Y_O\\\mathbf Y_M\end{pmatrix}
\sim\mathcal N\left[
\mu\begin{pmatrix}\mathbf1_O\\\mathbf1_M\end{pmatrix},
\sigma^2\begin{pmatrix}C_{OO}&C_{OM}\\C_{MO}&C_{MM}\end{pmatrix}
\right].
$$

Integrating out $\mathbf Y_M$ leaves

$$
\mathbf Y_O\sim\mathcal N(\mu\mathbf1_O,\sigma^2C_{OO}).
$$

Writing $n_O$ for the number observed, the ML estimates are therefore

$$
\widehat\mu=\frac{\mathbf1_O^TC_{OO}^{-1}\mathbf y_O}
{\mathbf1_O^TC_{OO}^{-1}\mathbf1_O},\qquad
\widehat\sigma^2=\frac{(\mathbf y_O-\widehat\mu\mathbf1_O)^TC_{OO}^{-1}
(\mathbf y_O-\widehat\mu\mathbf1_O)}{n_O}.
$$

Obuina constructs the observed-descendant tree while retaining the original root and all relevant path lengths, then delegates its Brownian fit/reconstruction to curupixa. This preserves the original $C_{OO}$, including shared paths above the observed tips' MRCA. Missing-only subtrees contribute no likelihood term. The calculation remains linear in the original node count; it does not construct the dense matrix shown above.

At least two distinct finite observed values, positive original edges and a positive fitted diffusion variance are required. An all-missing or one-observation trait cannot supply a finite positive-variance ML fit here. A large number of missing tips does not increase the variance estimator's denominator. The C result's `tips` field reports $n_O$, while `nodes` still refers to the full original tree.

## 3. Brownian prediction includes root-state uncertainty

For a fixed root mean $\mu$, Gaussian conditioning gives

$$
E(\mathbf Y_M\mid\mathbf y_O,\mu)=
\mu\mathbf1_M+C_{MO}C_{OO}^{-1}(\mathbf y_O-\mu\mathbf1_O),
$$

$$
\operatorname{Cov}(\mathbf Y_M\mid\mathbf y_O,\mu)=
\sigma^2(C_{MM}-C_{MO}C_{OO}^{-1}C_{OM}).
$$

For reconstruction, Obuina integrates the ancestral root state under a flat prior, with $\widehat\sigma^2$ plugged in. Define

$$
a=\mathbf1_O^TC_{OO}^{-1}\mathbf1_O,\qquad
\mathbf h=\mathbf1_M-C_{MO}C_{OO}^{-1}\mathbf1_O.
$$

The full missing-tip prediction is

$$
E(\mathbf Y_M\mid\mathbf y_O)=
\widehat\mu\mathbf1_M+C_{MO}C_{OO}^{-1}
(\mathbf y_O-\widehat\mu\mathbf1_O),
$$

$$
\operatorname{Cov}(\mathbf Y_M\mid\mathbf y_O)=\widehat\sigma^2
\left[C_{MM}-C_{MO}C_{OO}^{-1}C_{OM}+\frac{\mathbf h\mathbf h^T}{a}\right].
$$

The final term accounts for root-state uncertainty. It would be absent if the root were fixed at its fitted mean. Obuina returns each node's marginal variance rather than the full dense missing-tip covariance. Those missing predictions can be correlated even though only separate variances are returned.

For a whole unobserved child subtree, the Markov property provides a simple propagation rule. If parent $p$ has conditional mean $m_p$ and variance $v_p$, its unobserved child $c$ across edge $t_c$ has

$$
m_c=m_p,\qquad v_c=v_p+\widehat\sigma^2t_c.
$$

This explains why uncertainty grows down an unobserved lineage. An exactly observed tip still has zero variance. The prediction variances exclude diffusion-rate estimation uncertainty, tree uncertainty and measurement error. A conditional Gaussian interval uses $m\pm1.96\sqrt v$ under those assumptions.

## 4. PIC and regression use the observed data directly

For PIC, removing a missing-only subtree and contracting a path with one observed descendant preserves the contrast variance accumulated along that path. The original tree must satisfy the method's binary/positive-edge requirements. With $n_O$ retained observations, the result has $n_O-1$ contrasts. Contrast node IDs map to the corresponding branching nodes in the original tree.

For PGLS, let $I$ contain tips with finite response and all finite predictors. Obuina fits

$$
\mathbf y_I=X_I\boldsymbol\beta+\boldsymbol\epsilon_I,
\qquad\boldsymbol\epsilon_I\sim\mathcal N(0,\sigma^2C_{II}).
$$

It retains the covariance relative to the original root rather than arbitrarily rerooting the retained taxa. The result's `tips` field counts retained rows. The requirements $|I|>p$, full design rank and positive residual variance still apply. The CLI's observed/used flags distinguish a known response omitted because one of its predictors is missing.

No predicted response or predictor is substituted into this fit. Missing predictors require a model for their distribution and their association with the response before likelihood-based joint imputation can be justified. Complete-case fitting is the current explicit policy; it is not generally unbiased under every missingness mechanism.

## 5. Several columns are independent analyses

After building from the project root:

```sh
mkdir -p results/missing
build/obuina brownian --tree examples/tree.nwk --data examples/traits-missing.tsv --traits value,second --missing estimate --prefix results/missing/brownian
```

Each selected trait produces its own files, for example `brownian.trait1.nodes.tsv` and `brownian.trait2.nodes.tsv`. `brownian.traits.tsv` records `trait`, `prefix` and `status`; numeric suffixes keep arbitrary column names out of filesystem paths. Each column has its own observed set, root mean, variance and likelihood. The plural form requires a prefix to keep outputs separate; the ordinary `--trait` form remains available. CLI batches continue after a column fails, record `error` for that entry and return a nonzero exit code if any fail. They do not roll back earlier outputs. Files from a previous run of a failed column can still exist, so consult the current manifest's status rather than treating file existence as success.

For regression, `--predictors size,temperature` selects several predictor columns with the singular `--predictor size` retained for compatibility. UPGMA can use the selected primary numeric feature plus several predictors as extra distance features, all complete and appropriately scaled.

The C function `bt_brownian_fit_many` accepts a node-major `[node][trait]` matrix and returns independent fits. Unlike the CLI's continue-and-report behaviour, its owned result has an all-or-nothing contract: if a trait fails, no partial owned batch result is returned. The Python `brownian_many` interface accepts a mapping from column names to taxon/value mappings and returns a mapping of names to ordinary per-trait fit dictionaries; a failing trait raises an error without returning a partial batch.

These independent models do not estimate cross-trait covariance. A future multivariate Brownian model might specify

$$
\operatorname{Cov}(Y_{ia},Y_{jb})=C_{ij}\Sigma_{ab},
$$

where $\Sigma$ represents evolutionary covariance among trait dimensions. Such a model could borrow information across columns when predicting missing cells. Independent batches currently borrow information only within each trait column, even when the columns look correlated.

## 6. Discrete observation uncertainty is available through the API

For tip $i$ in latent state $s$, supply an observation-likelihood weight

$$
w_i(s)=p(\text{observed metadata}_i\mid Z_i=s).
$$

The C function `bt_discrete_asr_likelihoods` accepts node-major `[node][state]` weights; internal-node rows are ignored. Python exposes `Tree.discrete_likelihoods` with an explicit state alphabet. Each tip vector must be finite, nonnegative and have at least one positive entry.

| Observation | Weight vector for states A, B, C |
|---|---|
| Exact state B | `[0, 1, 0]` |
| Completely unknown | `[1, 1, 1]` |
| Known to be A or C | `[1, 0, 1]` |
| A noisy observation compatible with several states | A supplied vector such as `[0.8, 0.15, 0.05]` under an explicit observation model |

The vectors are likelihoods of observations given states. They need not sum to one. A classifier's posterior state probabilities cannot generally be used as observation likelihoods without accounting for the prior used by that classifier. If $q(s)=p(s\mid d)$ was computed under known prior $\pi_{\mathrm{classifier}}(s)>0$, Bayes' rule gives weights proportional to $q(s)/\pi_{\mathrm{classifier}}(s)$, subject to compatibility of that observation model with the analysis.

For Python likelihood-vector results, `observed` records whether a vector was supplied rather than omitted/`None`; it does not measure informativeness. An explicitly supplied `[1, 1, 1]` vector is flagged as supplied evidence even though it has the same marginal effect as an unknown `None` vector.

Multiplying all entries of one tip's vector by the same positive constant leaves its reconstructed marginals and fitted rate unchanged, while adding the logarithm of that constant to the total log likelihood. An all-ones vector is the canonical unknown observation. Merely normalising it to sum to one would add an arbitrary likelihood offset.

Likelihood pruning starts from $D_i(s)=w_i(s)$, rather than a one-hot exact state or all-ones unknown. The full-data upward pass and marginal normalisation are otherwise the same. This API supports soft evidence without a new state-transition model. The present TSV CLI exposes exact states and unknown markers; it does not yet accept columns of arbitrary likelihood vectors. Parsimony's public Obuina API still accepts exact/unknown states rather than weighting histories by these probabilities.

## 7. Reproduce missing and uncertain observations

The [example metadata](../examples/traits-missing.tsv) deliberately omit D's row and leave different cells unknown in A, B and C. Each continuous analysis has enough observed, distinct values to fit. Run:

```sh
build/obuina brownian --tree examples/tree.nwk --data examples/traits-missing.tsv --trait value --missing estimate --prefix results/missing/value
build/obuina discrete --tree examples/tree.nwk --data examples/traits-missing.tsv --trait class --missing estimate --prefix results/missing/class
build/obuina pic --tree examples/tree.nwk --data examples/traits-missing.tsv --trait value --missing estimate --prefix results/missing/pic
PYTHONPATH=python python3 examples/missing_and_uncertain.py
```

For these CLI results, `.nodes.tsv` includes `observed` and `used_for_fit` flags, and `.summary.tsv` records `trait`, `missing_policy`, `observed_tips`, `missing_tips` and `used_tips`. A prediction at a missing Brownian/discrete tip appears in the ordinary estimate/probability columns. Retaining those flags prevents a later reader from confusing reconstructed metadata with observed measurements.

For the `value` example, the Brownian fit uses only A and C. Rounded output is:

| Tip | Observed? | Conditional mean | Conditional variance |
|---|---|---:|---:|
| A | Yes | 1.200000 | 0 |
| B | No: explicit `NA` | 1.409302 | 0.260654 |
| C | Yes | 2.100000 | 0 |
| D | No: omitted row | 1.932558 | 0.193191 |

The two unobserved tips receive estimates with uncertainty while the fitted observation count remains two. This small example demonstrates the computation; inference about evolutionary rates from only two observations is weak.

## 8. Continuous measurement uncertainty: proposed next step

Continuous missingness is implemented; noisy continuous measurements require an additional observation model. A useful first extension would take measured value $y_i$ and known standard error $s_i$ for a latent taxon state $z_i$:

$$
y_i\mid z_i\sim\mathcal N(z_i,s_i^2),\qquad
\mathbf y\sim\mathcal N(X\boldsymbol\beta,V),\qquad
V=\sigma^2C+\operatorname{diag}(s_i^2).
$$

This model would combine evolutionary variation and independent measurement error. The additive measurement variances are in squared trait units and do not scale with the diffusion rate. Fitting $\sigma^2$ would therefore generally need a scalar likelihood optimisation instead of the present closed-form SSE/divisor calculation. Latent-tip reconstructions would have positive uncertainty even where a noisy measurement is available.

Potential inputs and their meanings are:

| User information | Required interpretation |
|---|---|
| Known measurement SD or variance | Observation error for that measurement, with units explicit |
| SE of a taxon mean | Uncertainty of the mean; under appropriate independent sampling, SE equals sample SD divided by $\sqrt n$ |
| Replicate measurements | A model distinguishing within-taxon sampling variation from evolution of the taxon mean |
| An interval or detection limit | An interval/censoring likelihood, or a clearly justified approximation |
| Correlated measurement errors | An observation covariance matrix rather than independent diagonal errors |
| Uncertain predictors | An errors-in-variables or joint latent-trait model; response-error weights alone are insufficient |

An SD among individuals is not automatically the SE of their mean. An interval's interpretation also matters: a range, a confidence interval and a detection limit imply different likelihoods. Input precision should follow the observation-generating process rather than being chosen merely to obtain narrower reconstructed intervals.

For a future C API, a size-tagged options object could carry an optional per-node observation-variance array, an explicit missing-value mask, root treatment and a flag controlling reconstruction allocation. An all-zeros variance array should reproduce exact-observation results. The existing API would remain a convenience wrapper with those defaults. Python could accept a separate taxon-to-variance or taxon-to-SE mapping and require an explicit convention. These are design options, not currently accepted parameters.

## 9. Missingness and uncertainty in the method itself

Missing-data likelihood assumes that the unobserved-state mechanism can be ignored under the chosen model, or that the analysis conditions on the observed pattern appropriately. Missing completely at random is a sufficient simple case. More general ignorability depends on which variables and missingness information are modelled. If extreme states are systematically harder to measure, prediction using only the evolutionary process can be biased. The missing-data switch does not fit a selection/detection model.

Observation uncertainty and uncertainty in evolutionary parameters are separate. The current Brownian/discrete predictions condition on fitted evolutionary rates and one tree. Brownian reconstruction includes ancestral root-state uncertainty under its flat prior, while excluding the uncertainty from estimating the diffusion variance. Future profile likelihood, parametric bootstrap or multiple-imputation procedures could propagate parameter and tree uncertainty. Multiple imputation would need draws from an appropriate joint conditional distribution, preserving correlation among missing tips and, in a multivariate model, among dimensions.

Use the returned variances and probabilities when describing uncertainty. Preserve observed/used flags and the original missing pattern, and report the process, tree and missingness assumptions alongside the estimates.
