---
tags: [phylogenetics, ASR, PCM, CTMC, parsimony, tutorial]
---
# Discrete traits: likelihood, marginal states and parsimony

[Tutorial contents](asr-pcm-tutorial.md) · [Previous: continuous traits](tutorial-continuous.md) · [Next: extensions](tutorial-extensions.md)

Obuina offers two approaches to discrete ASR. A continuous-time Markov chain, CTMC, assigns probabilities to evolutionary changes along edges. Parsimony assigns costs and finds minimum-cost reconstructions. Both use dynamic programming, but their outputs have different interpretations. The current likelihood model treats states as unordered; parsimony supports both unordered and ordered states.

## 1. The equal-rates continuous-time Markov chain

Let a trait have $K\geq2$ possible states. A generator matrix $Q$ specifies instantaneous transition rates. For Obuina' equal-rates model,

$$
Q_{ij}=\frac{r}{K-1}\quad(i\ne j),\qquad Q_{ii}=-r.
$$

The rows sum to zero. Every state has total exit rate $r$, so the expected waiting time before a change is $1/r$. Obuina' `rate` therefore measures expected state changes per branch unit. Some software instead calls each individual off-diagonal rate $q$ the ER rate; in that convention $r=(K-1)q$. They agree for two states.

Along an edge of length $t$, the transition matrix is $P(t)=\exp(Qt)$. For this model there is a closed form:

$$
P_{ij}(t)=\frac1K+\left(\mathbf1_{\{i=j\}}-\frac1K\right)
e^{-Krt/(K-1)}.
$$

At $t=0$, transitions are the identity. As $t$ increases, each row approaches the uniform stationary frequencies $\pi_i=1/K$. Endpoint equality does not mean that no changes occurred: the process may leave a state and return before the edge ends.

### Assumptions

The model assumes a fixed known state alphabet, equal rates among all pairs, the same generator across the tree, a uniform stationary root distribution and independent evolution on descendant edges conditional on their parent's state. Unknown tips allow every state. Branch lengths are finite and nonnegative; polytomies, unary nodes and unequal tip depths are accepted.

The alphabet must include biologically possible unobserved states when they are part of the hypothesis. Adding states changes both equilibrium frequencies and transition probabilities, even if those states never occur in the observations. Arbitrary category names can be relabelled; the underlying alphabet size and state correspondence must stay fixed.

Obuina currently implements this single ER likelihood model. Broader Mk interfaces, such as [phytools' discrete model family](https://github.com/liamrevell/phytools/blob/master/man/fitMk.Rd), also support unequal and ordered rates. Those require additional generators and parameter estimation here.

## 2. Calculate the tip likelihood by pruning

At node $v$, define the downward likelihood

$$
D_v(i)=p(\text{observations below }v\mid Z_v=i).
$$

For an observed tip in state $a$, $D_v(i)=\mathbf1_{\{i=a\}}$. For an unknown tip, $D_v(i)=1$ for every $i$. Unknown data therefore contribute no preference among states; they are not coded as an additional biological state.

For an internal node,

$$
D_v(i)=\prod_{c\in\operatorname{children}(v)}M_{c\to v}(i),
\qquad
M_{c\to v}(i)=\sum_jP_{ij}(t_c)D_c(j).
$$

Each sum integrates over a child's unknown state. The product combines independent descendant branches conditional on the parent state. At root $\rho$,

$$
L(r)=\sum_i\pi_iD_\rho(i).
$$

This sums over all ancestral assignments without enumerating $K$ to the power of the number of internal nodes. Curupixa performs pruning in log space and supplies the full reconstruction API. Obuina wraps those results and its rate-fitting/application policies; it no longer maintains a separate upward reconstruction kernel.

For a general generator, a dense transition application costs $O(K^2)$ per edge. ER has the special form

$$
M(i)=aD(i)+(1-a)\frac{\sum_jD(j)}K,
\qquad a=e^{-Krt/(K-1)},
$$

so each edge needs $O(K)$ work. A fixed-rate reconstruction takes $O(NK)$ time and storage apart from input validation/mapping overhead. Log-space calculations prevent products across many tips from underflowing.

## 3. Recover full marginal probabilities with a second pass

A non-root downward likelihood omits observations outside that node's subtree. Introduce an outside message

$$
U_v(i)=p(Z_v=i,\text{observations outside }v\text{'s subtree}).
$$

Start with $U_\rho(i)=\pi_i$. For child $c$ of $v$,

$$
U_c(j)=\sum_iU_v(i)P_{ij}(t_c)
\prod_{b\in\operatorname{children}(v),\ b\ne c}M_{b\to v}(i).
$$

The child receives information from its parent, ancestors and sibling subtrees. Its full marginal is

$$
p(Z_v=i\mid\mathbf y,r,\mathcal T)
=\frac{U_v(i)D_v(i)}{\sum_jU_v(j)D_v(j)}.
$$

The implementation rescales messages by state-independent constants for numerical stability; the final normalisation removes those constants. It handles zero-probability contributions explicitly, avoiding invalid subtraction of infinite log terms.

### Marginal, joint and history outputs

Obuina returns each node's marginal probability vector. Its CLI `state` column chooses the largest entry, breaking exact ties by the lowest state index. A **joint MAP** method would instead maximise the probability of a complete ancestral assignment. Independent marginal modes need not maximise that joint probability: individually probable states can form a combination with low joint probability. A **stochastic map** would additionally sample the events within edges, including unobserved reversals.

The [PastML guide](https://pastml.pasteur.fr/help) distinguishes marginal modes, joint assignments and MPPA summaries that retain selected sets of probable states. Obuina currently exposes full marginals and their modes, with joint reconstruction and mapping deferred. These conditional probabilities require a root prior but no MCMC calculation. With a fitted rate, they condition on that rate and omit its estimation uncertainty.

Names used by packages can differ: for example, [ape's `ace` documentation](https://search.r-project.org/CRAN/refmans/ape/html/ace.html) describes a subtree-only calculation under one of its settings. To compare results, check precisely which tips contribute to each reported state distribution, as well as the root prior and rate convention.

### A two-tip example with unequal lengths

Take `(A:1,B:3)R;`, observed A=0, B=1, and $r=0.25$. Write $a=e^{-2r}=e^{-0.5}$ and $b=e^{-6r}=e^{-1.5}$. The likelihood is

$$
L=\frac{1-ab}{4}\approx0.216166,
$$

and the root probability is

$$
p(R=0\mid A=0,B=1)=\frac{(1+a)(1-b)}{2(1-ab)}\approx0.721705.
$$

The shorter path to A retains more information about the root under this model. The log likelihood is approximately $-1.531708$.

```python
from obuina import Tree

with Tree("(A:1,B:3)R;") as tree:
  fit = tree.discrete({"A": 0, "B": 1}, states=[0, 1], rate=0.25)
  print(fit["log_likelihood"])
  print(fit["nodes"][0]["probabilities"])
```

Run Python snippets from the project root with `PYTHONPATH=python` set. If both edges had length zero and the observations still differed, the likelihood would be zero and the conditional marginals undefined. Obuina returns an error in that case. With all observations unknown, marginals are uniform and the likelihood contains no rate information.

## 4. Fix or estimate the rate

From the [tutorial setup](asr-pcm-tutorial.md#4-set-up-the-worked-examples), run

```sh
build/obuina discrete --tree examples/tree.nwk --data examples/traits.tsv --trait class --rate 1 --prefix results/tutorial/discrete-fixed
build/obuina discrete --tree examples/tree.nwk --data examples/traits.tsv --trait class --fit-rate --prefix results/tutorial/discrete-fitted
```

Without `--fit-rate`, the rate is fixed at the supplied value; its default is 1. With fitting, Obuina prepares one curupixa likelihood binding and searches

$$
r\in\left[\frac{10^{-8}}{\bar t},\frac{100}{\bar t}\right],
$$

where $\bar t$ is the mean non-root edge length, including zero-length edges. It first evaluates a 65-point logarithmic grid, then refines the best neighbouring interval by golden-section search. This bounded search can return boundary or plateau estimates and is not a proof of a global optimum for every possible surface. All-zero edges cannot identify a rate, so fitting returns an error.

The supplied four-tip example gives the following at fixed rate 1:

| Quantity | Value |
|---|---:|
| Log likelihood | -2.787311 |
| Root probability of `high` | 0.488507 |
| Root probability of `low` | 0.511493 |

The CLI's inferred alphabet is lexicographic: `state_0=high`, `state_1=low`. An explicit Python alphabet `["low", "high"]` reverses the probability-column order while preserving the labelled result. Here the root is weakly resolved; reducing the vector to `low` discards most of the uncertainty.

The fitted-rate example can reach a high-rate likelihood plateau, where ancestral marginals approach the uniform frequencies. Its precise returned rate is then weakly informative and can vary with floating-point arithmetic or search details. Inspect the score and probabilities and, when inference depends on the rate, evaluate a profile over fixed rates. The scalar summary currently has no dedicated boundary/identifiability flag.

## 5. Coding missing, ordered and unobserved states

Text classes are convenient for unordered traits. To retain a declared alphabet including unobserved states in the CLI, use integer observations and `--states K`, with codes `0..K-1`. Python accepts an explicit list of labels.

Unknown observations allow every state and can be inferred using the rest of the tree. The C `bt_discrete_asr_likelihoods` and Python `discrete_likelihoods` APIs now accept node/state observation-likelihood vectors: one-hot entries represent exact states, all ones represent unknowns and multiple positive entries represent ambiguity or an explicit noisy observation model. They represent $p(\text{observation}\mid\text{state})$, with arbitrary positive overall scale per tip. Classifier posteriors require attention to the classifier's prior. See [observation uncertainty](missing-data-and-uncertainty.md#6-discrete-observation-uncertainty-is-available-through-the-api) for the interpretation and scale effect on likelihoods. The TSV interface still accepts exact/unknown states only.

Omitted metadata rows require `--missing estimate`; explicit discrete unknown markers remain accepted without that flag. Python omitted keys require `allow_missing=True`, while an explicit unknown `None` remains supported. Unknown data provide no model of why observations are missing; informative missingness needs additional treatment.

Ordering must be specified through a transition model or cost matrix. Under current ER likelihood, `0` can change directly to `2` with exactly the same rate as `0` to `1`. Under ordered parsimony, that change costs two steps. An adjacent-state CTMC, where direct `0` to `2` instantaneous transitions are forbidden but two successive changes remain possible, is a future model described in the [extensions chapter](tutorial-extensions.md#4-ordered-discrete-likelihood).

## 6. Parsimony and the Sankoff recurrence

Let $c(i,j)$ be the cost assigned to an edge whose parent and child have states $i,j$. A complete reconstruction $\mathbf z$ has score

$$
S(\mathbf z)=\sum_{(u,v)\in E}c(z_u,z_v),
$$

subject to observed tip states. Obuina implements

$$
c_{\mathrm{unordered}}(i,j)=\mathbf1_{\{i\ne j\}},\qquad
c_{\mathrm{ordered}}(i,j)=|i-j|.
$$

For a node fixed to state $i$, define $F_v(i)$ as the minimum cost in its descendant subtree. At a known tip, the observed state has cost zero and all other states have infinite cost. At an unknown tip, all states cost zero. For an internal node,

$$
F_v(i)=\sum_{c\in\operatorname{children}(v)}\min_j\{c(i,j)+F_c(j)\}.
$$

The optimum score is $S_{\min}=\min_iF_\rho(i)$. This is the Sankoff dynamic programme. A general dense cost matrix requires $O(NK^2)$ work. For unordered costs, a global minimum of the child vector gives the transform in $O(K)$; for $|i-j|$, two directional minimum sweeps do the same. Obuina now delegates those numerical passes to curupixa's reconstruction API, retaining $O(NK)$ time/storage for both implemented cases.

### Recover all optimal states and one coherent assignment

An outside-cost pass computes $G_v(i)$, the minimum cost elsewhere in the tree conditional on $z_v=i$. A state belongs to the global optimal set at node $v$ exactly when

$$
F_v(i)+G_v(i)=S_{\min}.
$$

Different nodes' sets describe different possible global reconstructions. Choosing a state independently from each set can fail to produce an optimal complete assignment. Obuina also performs a consistent traceback: choose a minimum-cost root state, then, given parent state $i$, choose a child state minimising $c(i,j)+F_c(j)$. Ties select the lowest state index.

The resulting `state` column is one coherent optimum. `optimal_states` retains the ambiguity across all optima. Neither output assigns probabilities or counts optimal scenarios. The [phangorn ASR interface](https://klausvigo.github.io/phangorn/reference/ancestral.pml.html) similarly highlights the non-uniqueness of parsimony reconstructions; its available ambiguity conventions extend beyond those implemented here.

### A difference caused by ordering

For `(A,B)R;` with A=0, B=2 and declared alphabet `[0,1,2]`, unordered parsimony has score 1 with root set $\{0,2\}$. Ordered parsimony has score 2 with root set $\{0,1,2\}$. The intermediate root state 1 is optimal under the ordered cost although it is absent at the tips. The lowest-index traceback chooses root 0 in both cases.

Under unordered costs the score counts changed edges in the selected reconstruction. Under ordered costs it counts summed level differences: one edge can contribute several steps. It does not locate those steps in time or supply a probability distribution for their number.

## 7. Run parsimony and compare the objectives

```sh
build/obuina parsimony --tree examples/tree.nwk --data examples/traits.tsv --trait class --prefix results/tutorial/parsimony-unordered
build/obuina parsimony --tree examples/tree.nwk --data examples/traits.tsv --trait ordered --ordered --states 3 --prefix results/tutorial/parsimony-ordered
```

The supplied example has unordered class score **2** and ordered-level score **3**. These are analyses of different trait columns with different cost assumptions, so the score difference does not compare model fit to the same observation model. To isolate ordering, analyse the same integer trait once with each cost convention.

| Feature | ER likelihood ASR | Ordered/unordered parsimony |
|---|---|---|
| Edge information | Length and transition process | Topology and declared costs |
| Ancestral alternatives | Summed over using probabilities | Minimised over using costs |
| Ambiguity | Marginal probabilities | Possible states across global optima |
| Parameters | Fixed or fitted rate; specified uniform root prior | Cost rule and state alphabet |
| Unknown tips | All states allowed | All states allowed |
| Direction/history | Node marginals; historical direction needs rooting | One minimum-cost assignment plus state sets |

Parsimony supports missing lengths and even ignores finite negative lengths accepted by the general tree parser. Likelihood requires nonnegative lengths. Negative lengths have no evolutionary-exposure interpretation for the CTMC. Parsimony also accepts a one-state alphabet, with score zero for a constant observed trait. Both methods treat polytomies as **hard**: all children descend directly from that vertex, with no search over hidden binary resolutions.

Minimum-change and likelihood reconstructions can agree in some regimes, but finite branch lengths and probabilities of multiple changes can produce different results. Comparing the likelihood's log score to a parsimony cost has no direct statistical meaning. Choosing between them should follow the evolutionary/transition assumptions and the uncertainty representation needed for the question.
