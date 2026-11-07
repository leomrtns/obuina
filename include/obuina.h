/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef OBUINA_H
#define OBUINA_H
#include <stddef.h>
#include <stdbool.h>
#include <tree.h>
#ifdef __cplusplus
extern "C" {
#endif

#define BT_VERSION "0.2.0"
/* All functions return recoverable curupixa statuses. Constructors set *out=NULL
 * on failure. Results own their arrays; free with the matching function (NULL is
 * safe). Inputs are borrowed only during a call. Independent calls are reentrant.
 * Trait arrays are indexed by stable curupixa node ID; internal entries are ignored.
 * Callers supply arrays of at least crpx_tree_node_count(tree) elements. */
const char *bt_version(void);
const char *bt_status_string(crpx_status status);

typedef struct {
  size_t nodes, tips;
  double minimum_tip_depth, maximum_tip_depth;
  bool ultrametric;
} bt_tree_info;
/* Require named, unique tips and complete finite, nonnegative edge lengths.
 * Root stems are ignored. tolerance is relative to max(1, maximum tip depth). */
crpx_status bt_tree_inspect(const crpx_tree *tree, double tolerance, bt_tree_info *out);
/* Rebuild at the midpoint of the weighted tip diameter. All taxon distances are
 * preserved. Node IDs can change: remap traits after this operation. */
crpx_status bt_tree_midpoint_root(const crpx_tree *tree, crpx_tree **out);
/* Extend terminal edges to the maximum current tip depth. This minimises added
 * length with internal edges fixed and terminal edges only allowed to increase.
 * Returns a clone. This numerical transformation does not estimate dated ages. */
crpx_status bt_tree_extend_tips(const crpx_tree *tree, crpx_tree **out, double *total_added);
/* Euclidean distances over complete finite row-major [taxa][features] metadata,
 * followed by curupixa UPGMA. Scale features before calling if units differ.
 * A tree inferred from the response cannot independently establish its signal. */
crpx_status bt_tree_upgma(const char *const *names, size_t taxa, size_t features,
                          const double *metadata, crpx_tree **out);

typedef struct {
  size_t nodes, tips;
  double root_mean, sigma2, log_likelihood;
  double *mean, *variance;
} bt_brownian_result;
/* Single-trait Brownian ML with a fitted root mean and diffusion variance.
 * NAN denotes an unobserved tip; infinities are rejected. Strictly positive
 * edges, >=2 observed tips and nonconstant observed data are required. tips in
 * the result counts observed tips. Fit marginalises missing tips; mean/variance
 * include predictions at all original nodes, without treating predictions as data.
 * Marginals use a flat root prior and plug-in sigma2; parameter uncertainty is
 * excluded. Works on non-ultrametric trees and polytomies. */
crpx_status bt_brownian_fit(const crpx_tree *tree, const double *values, bt_brownian_result **out);
void bt_brownian_free(bt_brownian_result *result);

typedef struct {
  size_t traits;
  bt_brownian_result **results;
} bt_brownian_batch_result;
/* Independent column fits, not a joint correlated-trait model. Input layout is
 * node-major [node][trait], NAN may vary by column. Internal entries ignored.
 * All-or-nothing ownership: any failed column frees earlier results and leaves
 * *out=NULL. Use separate calls when per-column failures need retaining. */
crpx_status bt_brownian_fit_many(const crpx_tree *tree, const double *values, size_t traits,
                                 bt_brownian_batch_result **out);
void bt_brownian_batch_free(bt_brownian_batch_result *result);

typedef struct {
  size_t count;
  crpx_node_id *node;
  double *contrast, *variance;
} bt_pic_result;
/* Standardised independent contrasts; strictly bifurcating tree, positive edges.
 * NAN tips are omitted by contracting their induced unary paths; original tree
 * must be binary and >=2 observations remain. No imputation is used. Sign follows
 * child order; contrast node IDs refer to original branching vertices. */
crpx_status bt_pic(const crpx_tree *tree, const double *values, bt_pic_result **out);
void bt_pic_free(bt_pic_result *result);

typedef struct {
  size_t tips, parameters;
  double *coefficients, *standard_errors;
  double sigma2, log_likelihood;
} bt_pgls_result;
/* Brownian GLS, fixed covariance from the supplied root and branch lengths.
 * predictors is node-major [node][columns], excluding the automatically added
 * intercept. NAN response or predictor entries omit that tip (complete cases);
 * infinities are errors. Result tips counts retained rows. Full-rank predictors,
 * positive edges and retained tips>columns+1 are required.
 * sigma2 is ML SSE/tips; coefficient SEs use SSE/(tips-parameters).
 * Positive residual variance is required (exact fits have a boundary ML estimate).
 * Dense reference solver: O(tips^3) time and O(tips^2) memory. */
crpx_status bt_pgls(const crpx_tree *tree, const double *response, const double *predictors,
                     size_t columns, bt_pgls_result **out);
void bt_pgls_free(bt_pgls_result *result);

typedef struct {
  size_t nodes, states;
  double rate, log_likelihood;
  double *probabilities; /* node-major [nodes][states], marginal probabilities */
} bt_discrete_result;
/* Equal-rates CTMC, uniform stationary root prior, rate in expected changes per
 * branch-length unit (curupixa normalisation). observations[node] is a state
 * index, SIZE_MAX denotes unknown. >=2 states; finite rate>0. If fit_rate=true,
 * bounded scalar ML searches [1e-8/mean_edge, 100/mean_edge]. Boundary estimates
 * are possible and should be inspected. No ordered likelihood model is implied. */
crpx_status bt_discrete_asr(const crpx_tree *tree, const size_t *observations, size_t states,
                            double rate, bool fit_rate, bt_discrete_result **out);
void bt_discrete_free(bt_discrete_result *result);
/* Uncertain discrete observations: node-major [node][state] likelihood weights,
 * proportional to p(observation | latent state), not necessarily state posteriors.
 * Tip entries finite/nonnegative with at least one positive; internal rows ignored.
 * All ones is unknown, one-hot is exact. Scaling a tip row shifts logL by log(scale)
 * and preserves marginals/fitted rate. The same rate/root conventions apply. */
crpx_status bt_discrete_asr_likelihoods(const crpx_tree *tree, const double *tip_likelihoods,
                                        size_t states, double rate, bool fit_rate,
                                        bt_discrete_result **out);

typedef struct {
  size_t nodes, states;
  double score;
  size_t *assignment; /* one coherent optimal reconstruction, ties choose lowest */
  unsigned char *optimal; /* [nodes][states]: occurs in some global optimum */
} bt_parsimony_result;
/* Hard-polytomy Sankoff reconstruction. Unordered cost is 0/1; ordered cost is
 * abs(i-j), so caller state coding defines distances. Branch lengths are ignored.
 * Unknown tips use SIZE_MAX. optimal is a set indicator, never a probability. */
crpx_status bt_parsimony_asr(const crpx_tree *tree, const size_t *observations, size_t states,
                             bool ordered, bt_parsimony_result **out);
void bt_parsimony_free(bt_parsimony_result *result);

#ifdef __cplusplus
}
#endif
#endif
