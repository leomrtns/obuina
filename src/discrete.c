/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "obuina.h"
#include <character_data.h>
#include <likelihood.h>
#include <reconstruction.h>
#include <substitution_model.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  crpx_tree *tree;
  crpx_character_data *data;
  crpx_substitution_model *model;
  crpx_likelihood_binding *binding;
  crpx_likelihood_workspace *likelihood;
  crpx_reconstruction_workspace *reconstruction;
  size_t nodes, states;
  double *lengths, *edges;
} discrete_workspace;

static void workspace_clear(discrete_workspace *w)
{
  crpx_likelihood_binding_free(w->binding);
  crpx_likelihood_workspace_free(w->likelihood);
  crpx_reconstruction_workspace_free(w->reconstruction);
  crpx_substitution_model_free(w->model);
  crpx_character_data_free(w->data);
  crpx_tree_free(w->tree);
  free(w->lengths); free(w->edges);
}

/* Synthetic labels on a private clone adapt node-indexed arrays to curupixa's
 * named-character interface without altering labels, topology or stable IDs. */
static crpx_status prepare(const crpx_tree *tree, const size_t *observations,
                            const double *tip_likelihoods, size_t k, bool lengths,
                            discrete_workspace *w)
{
  memset(w, 0, sizeof(*w));
  if (!tree || (!observations && !tip_likelihoods) || !k || (lengths && k < 2))
    return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_tree_validate(tree);
  if (status != CRPX_OK) return status;
  w->nodes = crpx_tree_node_count(tree); w->states = k;
  size_t n = w->nodes, tips = crpx_tree_leaf_count(tree);
  if (!n || !tips) return CRPX_ERROR_INVALID_ARGUMENT;
  if (k > SIZE_MAX / n || n * k > SIZE_MAX / sizeof(double) ||
      tips > SIZE_MAX / sizeof(char *) || tips > SIZE_MAX / sizeof(size_t)) return CRPX_ERROR_OVERFLOW;
  status = crpx_tree_clone(tree, &w->tree);
  if (status != CRPX_OK) return status;
  const char **names = calloc(tips, sizeof(*names));
  size_t *states = observations ? malloc(tips * sizeof(*states)) : NULL;
  double *values = tip_likelihoods ? malloc(tips * k * sizeof(*values)) : NULL;
  if (!names || (observations && !states) || (tip_likelihoods && !values)) {
    status = CRPX_ERROR_ALLOCATION; goto done;
  }
  if (lengths) {
    w->lengths = calloc(n, sizeof(double)); w->edges = calloc(n, sizeof(double));
    if (!w->lengths || !w->edges) { status = CRPX_ERROR_ALLOCATION; goto done; }
  }
  size_t tip = 0;
  for (size_t i = 0; i < n; ++i) {
    if (lengths && (crpx_node_id)i != crpx_tree_root(tree)) {
      status = crpx_tree_branch_length(tree, (crpx_node_id)i, w->lengths + i);
      if (status != CRPX_OK) goto done;
      if (!isfinite(w->lengths[i]) || w->lengths[i] < 0) {
        status = CRPX_ERROR_INVALID_ARGUMENT; goto done;
      }
    }
    if (crpx_tree_first_child(tree, (crpx_node_id)i) != CRPX_NODE_NONE) continue;
    char name[32];
    snprintf(name, sizeof(name), "%zu", i);
    status = crpx_tree_set_label(w->tree, (crpx_node_id)i, name);
    if (status != CRPX_OK) goto done;
    names[tip] = crpx_tree_label(w->tree, (crpx_node_id)i);
    if (observations) {
      if (observations[i] != SIZE_MAX && observations[i] >= k) {
        status = CRPX_ERROR_OUT_OF_RANGE; goto done;
      }
      states[tip] = observations[i];
    } else memcpy(values + tip * k, tip_likelihoods + i * k, k * sizeof(double));
    ++tip;
  }
  if (observations)
    status = crpx_character_data_from_states(names, tips, 1, k, states, NULL, false, &w->data);
  else
    status = crpx_character_data_from_tip_likelihoods(names, tips, 1, k, values, NULL, false, &w->data);
  if (status == CRPX_OK) status = crpx_reconstruction_workspace_new(&w->reconstruction);
  if (status == CRPX_OK && lengths) status = crpx_substitution_model_poisson(k, &w->model);
done:
  free(names); free(states); free(values);
  return status;
}

static void set_rate(discrete_workspace *w, double rate, bool retain, crpx_likelihood_options *options)
{
  for (size_t i = 0; i < w->nodes; ++i) {
    /* Saturate only beyond the representable ER long-edge limit. */
    long double length = (long double)rate * w->lengths[i];
    w->edges[i] = length > DBL_MAX ? DBL_MAX : (double)length;
  }
  crpx_likelihood_options_init(options);
  options->edge_values = w->edges; options->edge_count = w->nodes;
  options->retain_partials = retain;
}

static crpx_status evaluate(discrete_workspace *w, double rate, double *log_likelihood)
{
  crpx_likelihood_options options;
  set_rate(w, rate, false, &options);
  return crpx_likelihood_evaluate_prepared(w->binding, w->model, &options, w->likelihood, log_likelihood);
}

/* One prepared binding supplies every candidate likelihood. A logarithmic grid
 * identifies a basin; golden-section refinement improves it, including bounds. */
static crpx_status fit_rate(discrete_workspace *w, double *rate)
{
  long double total = 0;
  for (size_t i = 0; i < w->nodes; ++i) total += w->lengths[i];
  if (total <= 0 || w->nodes < 2) return CRPX_ERROR_INVALID_ARGUMENT;
  double mean = (double)(total / (w->nodes - 1));
  if (!(mean > 0) || !isfinite(mean)) return CRPX_ERROR_OVERFLOW;
  double lower = log(1e-8) - log(mean), upper = log(100.0) - log(mean);
  if (lower < log(DBL_TRUE_MIN) || upper > log(DBL_MAX)) return CRPX_ERROR_OVERFLOW;
  const size_t intervals = 64;
  size_t best = 0;
  double best_log = -INFINITY, best_x = lower;
  for (size_t i = 0; i <= intervals; ++i) {
    double x = lower + (upper - lower) * (double)i / intervals, ll;
    crpx_status status = evaluate(w, exp(x), &ll);
    if (status != CRPX_OK) return status;
    if (ll > best_log) { best_log = ll; best_x = x; best = i; }
  }
  if (best_log == -INFINITY) return CRPX_ERROR_INVALID_ARGUMENT;
  double a = lower + (upper - lower) * (double)(best ? best - 1 : 0) / intervals;
  double b = lower + (upper - lower) * (double)(best < intervals ? best + 1 : intervals) / intervals;
  const double ratio = 0.6180339887498948482;
  double x1 = b - ratio * (b - a), x2 = a + ratio * (b - a), f1, f2;
  crpx_status status = evaluate(w, exp(x1), &f1);
  if (status == CRPX_OK) status = evaluate(w, exp(x2), &f2);
  if (status != CRPX_OK) return status;
  for (size_t iteration = 0; iteration < 80 && b - a > 1e-10; ++iteration) {
    if (f1 > best_log) { best_log = f1; best_x = x1; }
    if (f2 > best_log) { best_log = f2; best_x = x2; }
    if (f1 >= f2) {
      b = x2; x2 = x1; f2 = f1; x1 = b - ratio * (b - a);
      status = evaluate(w, exp(x1), &f1);
    } else {
      a = x1; x1 = x2; f1 = f2; x2 = a + ratio * (b - a);
      status = evaluate(w, exp(x2), &f2);
    }
    if (status != CRPX_OK) return status;
  }
  *rate = exp(best_x);
  return CRPX_OK;
}

void bt_discrete_free(bt_discrete_result *result)
{ if (result) { free(result->probabilities); free(result); } }

static crpx_status reconstruct(discrete_workspace *w, double rate, bool estimate_rate, bt_discrete_result **out)
{
  crpx_status status = CRPX_OK;
  bt_discrete_result *result = NULL;
  crpx_likelihood_reconstruction *upstream = NULL;
  if (estimate_rate) {
    status = crpx_likelihood_workspace_new(&w->likelihood);
    if (status == CRPX_OK) status = crpx_likelihood_binding_new(w->tree, w->data, &w->binding);
    if (status == CRPX_OK) status = fit_rate(w, &rate);
    if (status != CRPX_OK) goto done;
  }
  crpx_likelihood_options options;
  set_rate(w, rate, true, &options);
  status = crpx_likelihood_reconstruct(w->tree, w->data, w->model, &options, w->reconstruction, &upstream);
  if (status != CRPX_OK) goto done;
  result = calloc(1, sizeof(*result));
  if (!result) { status = CRPX_ERROR_ALLOCATION; goto done; }
  result->nodes = w->nodes; result->states = w->states; result->rate = rate;
  result->log_likelihood = crpx_likelihood_reconstruction_score(upstream);
  result->probabilities = malloc(w->nodes * w->states * sizeof(double));
  if (!result->probabilities) { status = CRPX_ERROR_ALLOCATION; goto done; }
  for (size_t node = 0; node < w->nodes; ++node)
    memcpy(result->probabilities + node * w->states,
      crpx_likelihood_reconstruction_probabilities(upstream, (crpx_node_id)node, 0),
      w->states * sizeof(double));
  *out = result; result = NULL;
done:
  crpx_likelihood_reconstruction_free(upstream); bt_discrete_free(result);
  return status;
}

crpx_status bt_discrete_asr(const crpx_tree *tree, const size_t *observations, size_t k,
                            double rate, bool estimate_rate, bt_discrete_result **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!isfinite(rate) || rate <= 0) return CRPX_ERROR_INVALID_ARGUMENT;
  discrete_workspace w;
  crpx_status status = prepare(tree, observations, NULL, k, true, &w);
  if (status == CRPX_OK) status = reconstruct(&w, rate, estimate_rate, out);
  workspace_clear(&w);
  return status;
}

crpx_status bt_discrete_asr_likelihoods(const crpx_tree *tree, const double *tip_likelihoods, size_t k,
                                       double rate, bool estimate_rate, bt_discrete_result **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!isfinite(rate) || rate <= 0) return CRPX_ERROR_INVALID_ARGUMENT;
  discrete_workspace w;
  crpx_status status = prepare(tree, NULL, tip_likelihoods, k, true, &w);
  if (status == CRPX_OK) status = reconstruct(&w, rate, estimate_rate, out);
  workspace_clear(&w);
  return status;
}

void bt_parsimony_free(bt_parsimony_result *result)
{ if (result) { free(result->assignment); free(result->optimal); free(result); } }

crpx_status bt_parsimony_asr(const crpx_tree *tree, const size_t *observations, size_t k,
                             bool ordered, bt_parsimony_result **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  discrete_workspace w;
  crpx_status status = prepare(tree, observations, NULL, k, false, &w);
  bt_parsimony_result *result = NULL;
  crpx_parsimony_reconstruction *upstream = NULL;
  if (status != CRPX_OK) goto done;
  if ((long double)(w.nodes - 1) * (ordered ? k - 1 : 1) > 9007199254740991.0L) {
    status = CRPX_ERROR_OVERFLOW; goto done;
  }
  status = crpx_parsimony_reconstruct(w.tree, w.data,
    ordered ? CRPX_PARSIMONY_ORDERED : CRPX_PARSIMONY_UNORDERED, NULL, 0, 0, w.reconstruction, &upstream);
  if (status != CRPX_OK) goto done;
  result = calloc(1, sizeof(*result));
  if (!result) { status = CRPX_ERROR_ALLOCATION; goto done; }
  result->nodes = w.nodes; result->states = k;
  result->score = crpx_parsimony_reconstruction_score(upstream);
  result->assignment = malloc(w.nodes * sizeof(size_t));
  result->optimal = malloc(w.nodes * k * sizeof(unsigned char));
  if (!result->assignment || !result->optimal) { status = CRPX_ERROR_ALLOCATION; goto done; }
  for (size_t node = 0; node < w.nodes; ++node) {
    result->assignment[node] = crpx_parsimony_reconstruction_assignment(upstream, (crpx_node_id)node, 0);
    memcpy(result->optimal + node * k,
      crpx_parsimony_reconstruction_optimal_states(upstream, (crpx_node_id)node, 0), k);
  }
  *out = result; result = NULL;
done:
  crpx_parsimony_reconstruction_free(upstream); bt_parsimony_free(result); workspace_clear(&w);
  return status;
}
