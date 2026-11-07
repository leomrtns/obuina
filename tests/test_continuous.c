/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "obuina.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void close_to(double actual, double expected)
{
  if (!(fabs(actual - expected) <= 2e-10 * fmax(1.0, fabs(expected)))) {
    fprintf(stderr, "Expected %.17g, received %.17g\n", expected, actual);
    abort();
  }
}

static crpx_tree *parse(const char *newick)
{
  crpx_tree *tree = NULL;
  size_t offset;
  assert(crpx_tree_from_newick(newick, &tree, &offset) == CRPX_OK);
  return tree;
}

static crpx_node_id find(const crpx_tree *tree, const char *label)
{
  for (size_t i = 0; i < crpx_tree_node_count(tree); ++i) {
    const char *name = crpx_tree_label(tree, (crpx_node_id)i);
    if (name && !strcmp(name, label)) return (crpx_node_id)i;
  }
  abort();
}

static double *observations(const crpx_tree *tree, const char *const *names, const double *values, size_t n)
{
  double *data = malloc(crpx_tree_node_count(tree) * sizeof(double));
  assert(data);
  for (size_t i = 0; i < crpx_tree_node_count(tree); ++i) data[i] = NAN;
  for (size_t i = 0; i < n; ++i) data[find(tree, names[i])] = values[i];
  return data;
}

static void brownian_dense_reference(void)
{
  crpx_tree *tree = parse("((A:1,B:2)I:1,C:4)R;");
  const char *names[] = {"A", "B", "C"};
  double values[] = {1, 4, 7};
  double *data = observations(tree, names, values, 3);
  bt_brownian_result *result = NULL;
  assert(bt_brownian_fit(tree, data, &result) == CRPX_OK);
  /* Independently calculated C=[[2,1,0],[1,3,0],[0,0,4]], det(C)=20,
   * C^-1=block([[3,-1],[-1,2]]/5,1/4). */
  double precision = 0.4 + 0.2 + 0.25;
  double root = (0.4 + 0.8 + 1.75) / precision;
  double sse = 17.65 - 2.95 * 2.95 / precision;
  double sigma2 = sse / 3.0;
  close_to(result->root_mean, root);
  close_to(result->sigma2, sigma2);
  close_to(result->log_likelihood, -0.5 * (3.0 * (log(2.0 * acos(-1.0) * sigma2) + 1.0) + log(20.0)));
  close_to(result->variance[crpx_tree_root(tree)], sigma2 / precision);
  /* With flat root prior, root and I have joint precision [[1.25,-1],[-1,2.5]]
   * and information [7/4,1+4/2]. Its inverse gives both node marginals. */
  double determinant = 1.25 * 2.5 - 1.0;
  close_to(result->mean[find(tree, "I")], (1.75 + 1.25 * 3.0) / determinant);
  close_to(result->variance[find(tree, "I")], sigma2 * 1.25 / determinant);
  for (size_t i = 0; i < 3; ++i) {
    close_to(result->mean[find(tree, names[i])], values[i]);
    close_to(result->variance[find(tree, names[i])], 0.0);
  }
  bt_pgls_result *gls = NULL;
  assert(bt_pgls(tree, data, NULL, 0, &gls) == CRPX_OK);
  close_to(gls->coefficients[0], root);
  close_to(gls->sigma2, result->sigma2);
  close_to(gls->log_likelihood, result->log_likelihood);
  close_to(gls->standard_errors[0], sqrt(sse / 2.0 / precision));
  bt_pgls_free(gls); bt_brownian_free(result);
  free(data); crpx_tree_free(tree);
}

static void polytomy_and_unary(void)
{
  crpx_tree *tree = parse("((A:1,B:2,C:4)I:0.75)R;");
  const char *names[] = {"A", "B", "C"};
  double values[] = {1, 4, 7};
  double *data = observations(tree, names, values, 3);
  bt_brownian_result *result = NULL;
  assert(bt_brownian_fit(tree, data, &result) == CRPX_OK);
  double mean = (1.0 + 4.0 / 2.0 + 7.0 / 4.0) / 1.75;
  double sse = (1.0 - mean) * (1.0 - mean) +
    (4.0 - mean) * (4.0 - mean) / 2.0 + (7.0 - mean) * (7.0 - mean) / 4.0;
  close_to(result->root_mean, mean);
  close_to(result->sigma2, sse / 3.0);
  close_to(result->mean[find(tree, "I")], mean);
  close_to(result->variance[find(tree, "I")], result->sigma2 / 1.75);
  close_to(result->variance[find(tree, "R")], result->sigma2 * (1.0 / 1.75 + 0.75));
  bt_pgls_result *gls = NULL;
  assert(bt_pgls(tree, data, NULL, 0, &gls) == CRPX_OK);
  close_to(result->log_likelihood, gls->log_likelihood);
  bt_pgls_free(gls); bt_brownian_free(result); free(data); crpx_tree_free(tree);
}

static void contrasts(void)
{
  crpx_tree *tree = parse("((A:1,B:1)I:1,C:2)R;");
  const char *names[] = {"A", "B", "C"};
  double values[] = {1, 3, 7};
  double *data = observations(tree, names, values, 3);
  bt_pic_result *result = NULL;
  assert(bt_pic(tree, data, &result) == CRPX_OK);
  assert(result->count == 2);
  assert(result->node[0] == find(tree, "I"));
  close_to(result->contrast[0], -2.0 / sqrt(2.0));
  close_to(result->variance[0], 2.0);
  close_to(result->contrast[1], -5.0 / sqrt(3.5));
  close_to(result->variance[1], 3.5);
  bt_pic_free(result); free(data); crpx_tree_free(tree);
  tree = parse("(A:1,B:1,C:2);");
  data = observations(tree, names, values, 3);
  result = (bt_pic_result *)1;
  assert(bt_pic(tree, data, &result) == CRPX_ERROR_TOPOLOGY && !result);
  free(data); crpx_tree_free(tree);
}

static void regression(void)
{
  crpx_tree *tree = parse("(A:1,B:1,C:1,D:1);");
  const char *names[] = {"A", "B", "C", "D"};
  double values[] = {1, 2, 4, 8}, x_values[] = {0, 1, 2, 3};
  double *data = observations(tree, names, values, 4);
  double *x = observations(tree, names, x_values, 4);
  bt_pgls_result *result = NULL;
  assert(bt_pgls(tree, data, x, 1, &result) == CRPX_OK);
  assert(result->tips == 4 && result->parameters == 2);
  close_to(result->coefficients[0], 0.3);
  close_to(result->coefficients[1], 2.3);
  close_to(result->sigma2, 2.3 / 4.0);
  close_to(result->standard_errors[0], sqrt(0.805));
  close_to(result->standard_errors[1], sqrt(0.23));
  bt_pgls_free(result);
  /* Unit changes in predictor coding preserve fitted values and uncertainty. */
  for (size_t i = 0; i < crpx_tree_node_count(tree); ++i) x[i] *= 1000.0;
  assert(bt_pgls(tree, data, x, 1, &result) == CRPX_OK);
  close_to(result->coefficients[0], 0.3);
  close_to(result->coefficients[1], 0.0023);
  close_to(result->standard_errors[1], sqrt(0.23) / 1000.0);
  bt_pgls_free(result);
  for (size_t i = 0; i < crpx_tree_node_count(tree); ++i) x[i] = 1.0;
  result = (bt_pgls_result *)1;
  assert(bt_pgls(tree, data, x, 1, &result) == CRPX_ERROR_INVALID_ARGUMENT && !result);
  x[find(tree, "A")] = NAN;
  assert(bt_pgls(tree, data, x, 1, &result) == CRPX_ERROR_INVALID_ARGUMENT && !result);
  free(x); free(data); crpx_tree_free(tree);
}


/* Independent dense Gaussian precision calculation for all unobserved states,
 * including missing tips: every edge contributes its quadratic penalty. */
static void dense_internal_marginals(const crpx_tree *tree, const double *values,
                                     const bt_brownian_result *result)
{
  size_t internal_count = 0;
  size_t index[32];
  size_t nodes = crpx_tree_node_count(tree);
  assert(nodes <= 32);
  for (size_t i = 0; i < nodes; ++i) {
    bool observed = crpx_tree_first_child(tree, (crpx_node_id)i) == CRPX_NODE_NONE && isfinite(values[i]);
    index[i] = observed ? SIZE_MAX : internal_count++;
  }
  assert(internal_count <= 16);
  double augmented[16][32] = {{0}}, information[16] = {0};
  for (size_t node = 0; node < nodes; ++node) {
    crpx_node_id parent = crpx_tree_parent(tree, (crpx_node_id)node);
    if (parent == CRPX_NODE_NONE) continue;
    double length;
    assert(crpx_tree_branch_length(tree, (crpx_node_id)node, &length) == CRPX_OK);
    double precision = 1.0 / length;
    size_t p = index[parent], c = index[node];
    augmented[p][p] += precision;
    if (c == SIZE_MAX) information[p] += precision * values[node];
    else {
      augmented[c][c] += precision;
      augmented[p][c] -= precision;
      augmented[c][p] -= precision;
    }
  }
  for (size_t i = 0; i < internal_count; ++i) augmented[i][internal_count + i] = 1.0;
  for (size_t i = 0; i < internal_count; ++i) {
    double diagonal = augmented[i][i];
    assert(diagonal > 0.0);
    for (size_t k = 0; k < 2 * internal_count; ++k) augmented[i][k] /= diagonal;
    for (size_t j = 0; j < internal_count; ++j) {
      if (j == i) continue;
      double multiplier = augmented[j][i];
      for (size_t k = 0; k < 2 * internal_count; ++k) augmented[j][k] -= multiplier * augmented[i][k];
    }
  }
  for (size_t node = 0; node < nodes; ++node) {
    if (index[node] == SIZE_MAX) continue;
    size_t i = index[node];
    double mean = 0.0;
    for (size_t j = 0; j < internal_count; ++j) mean += augmented[i][internal_count + j] * information[j];
    close_to(result->mean[node], mean);
    close_to(result->variance[node], result->sigma2 * augmented[i][internal_count + i]);
  }
}

static void nested_marginals(void)
{
  crpx_tree *tree = parse("(((A:0.2,B:3,C:0.7):2):0.3,(D:0.1,E:0.8):1.5,F:6);");
  const char *names[] = {"A", "B", "C", "D", "E", "F"};
  double values[] = {2, -1, 5, 3, 8, -2};
  double *data = observations(tree, names, values, 6);
  bt_brownian_result *result = NULL;
  assert(bt_brownian_fit(tree, data, &result) == CRPX_OK);
  dense_internal_marginals(tree, data, result);
  bt_pgls_result *gls = NULL;
  assert(bt_pgls(tree, data, NULL, 0, &gls) == CRPX_OK);
  close_to(gls->coefficients[0], result->root_mean);
  close_to(gls->sigma2, result->sigma2);
  close_to(gls->log_likelihood, result->log_likelihood);
  bt_pgls_free(gls); bt_brownian_free(result); free(data); crpx_tree_free(tree);
}


static void multiple_predictors(void)
{
  crpx_tree *tree = parse("(A:1,B:1,C:1,D:1,E:1);");
  const char *names[] = {"A", "B", "C", "D", "E"};
  double values[] = {-1.9, 1.1, 4.6, 5.1, 6.1};
  double first[] = {-2, -1, 0, 1, 2}, second[] = {2, -1, -2, -1, 2};
  double *data = observations(tree, names, values, 5);
  double *predictors = calloc(crpx_tree_node_count(tree) * 2, sizeof(double));
  assert(predictors);
  for (size_t i = 0; i < 5; ++i) {
    size_t node = find(tree, names[i]);
    predictors[node * 2] = first[i];
    predictors[node * 2 + 1] = second[i];
  }
  bt_pgls_result *result = NULL;
  assert(bt_pgls(tree, data, predictors, 2, &result) == CRPX_OK);
  close_to(result->coefficients[0], 3.0);
  close_to(result->coefficients[1], 2.0);
  close_to(result->coefficients[2], -0.5);
  close_to(result->sigma2, 0.14);
  close_to(result->standard_errors[0], sqrt(0.07));
  close_to(result->standard_errors[1], sqrt(0.035));
  close_to(result->standard_errors[2], sqrt(0.025));
  bt_pgls_free(result); free(data); free(predictors); crpx_tree_free(tree);
}

static bool ancestor(const crpx_tree *tree, crpx_node_id node, crpx_node_id candidate)
{
  for (; node != CRPX_NODE_NONE; node = crpx_tree_parent(tree, node)) if (node == candidate) return true;
  return false;
}

/* Dense observed-tip covariance is independent of subtree construction. */
static void dense_observed_likelihood(const crpx_tree *tree, const double *values,
                                       const bt_brownian_result *fit)
{
  crpx_node_id tips[16];
  size_t count = 0, nodes = crpx_tree_node_count(tree);
  for (size_t i = 0; i < nodes; ++i)
    if (crpx_tree_first_child(tree, (crpx_node_id)i) == CRPX_NODE_NONE && isfinite(values[i]))
      tips[count++] = (crpx_node_id)i;
  assert(count <= 16 && fit->tips == count);
  double matrix[16][32] = {{0}}, logdet = 0;
  for (size_t i = 0; i < count; ++i) {
    for (size_t j = 0; j < count; ++j) {
      for (size_t node = 0; node < nodes; ++node) {
        if ((crpx_node_id)node == crpx_tree_root(tree)) continue;
        if (ancestor(tree, tips[i], (crpx_node_id)node) && ancestor(tree, tips[j], (crpx_node_id)node)) {
          double length;
          assert(crpx_tree_branch_length(tree, (crpx_node_id)node, &length) == CRPX_OK);
          matrix[i][j] += length;
        }
      }
    }
    matrix[i][count + i] = 1;
  }
  for (size_t i = 0; i < count; ++i) {
    double diagonal = matrix[i][i];
    assert(diagonal > 0); logdet += log(diagonal);
    for (size_t k = 0; k < 2 * count; ++k) matrix[i][k] /= diagonal;
    for (size_t j = 0; j < count; ++j) if (j != i) {
      double multiplier = matrix[j][i];
      for (size_t k = 0; k < 2 * count; ++k) matrix[j][k] -= multiplier * matrix[i][k];
    }
  }
  double precision = 0, information = 0, quadratic = 0;
  for (size_t i = 0; i < count; ++i)
    for (size_t j = 0; j < count; ++j) {
      double inverse = matrix[i][count + j];
      precision += inverse;
      information += inverse * values[tips[j]];
      quadratic += values[tips[i]] * inverse * values[tips[j]];
    }
  double mean = information / precision;
  double sigma2 = (quadratic - information * mean) / count;
  close_to(fit->root_mean, mean); close_to(fit->sigma2, sigma2);
  close_to(fit->log_likelihood, -0.5 * (count * (log(2 * acos(-1.0) * sigma2) + 1) + logdet));
}

static void missing_brownian(void)
{
  crpx_tree *tree = parse("(((A:1,B:2)I:1,C:3)J:2,(D:1,E:2)K:4)R;");
  const char *names[] = {"A", "B", "C"};
  double values[] = {1, 4, 7};
  double *data = observations(tree, names, values, 3);
  bt_brownian_result *result = NULL;
  assert(bt_brownian_fit(tree, data, &result) == CRPX_OK);
  assert(result->nodes == crpx_tree_node_count(tree) && result->tips == 3);
  dense_observed_likelihood(tree, data, result);
  dense_internal_marginals(tree, data, result);
  crpx_node_id root = find(tree, "R"), missing = find(tree, "D");
  close_to(result->mean[missing], result->mean[root]);
  close_to(result->variance[missing], result->variance[root] + result->sigma2 * 5);
  bt_brownian_free(result);
  data[find(tree, "B")] = NAN;
  assert(bt_brownian_fit(tree, data, &result) == CRPX_OK);
  assert(result->tips == 2);
  dense_observed_likelihood(tree, data, result); dense_internal_marginals(tree, data, result);
  assert(result->variance[find(tree, "B")] > 0);
  bt_brownian_free(result);
  data[find(tree, "A")] = NAN;
  assert(bt_brownian_fit(tree, data, &result) == CRPX_ERROR_INVALID_ARGUMENT && !result);
  data[find(tree, "C")] = NAN;
  assert(bt_brownian_fit(tree, data, &result) == CRPX_ERROR_INVALID_ARGUMENT && !result);
  data[find(tree, "A")] = 1; data[find(tree, "C")] = 1;
  assert(bt_brownian_fit(tree, data, &result) == CRPX_ERROR_INVALID_ARGUMENT && !result);
  data[find(tree, "A")] = INFINITY;
  assert(bt_brownian_fit(tree, data, &result) == CRPX_ERROR_INVALID_ARGUMENT && !result);
  free(data); crpx_tree_free(tree);
}

static void missing_contrasts(void)
{
  crpx_tree *tree = parse("(((A:1,B:2)I:3,C:4)J:2,(D:1,E:2)K:5)R;");
  const char *names[] = {"A", "C", "E"};
  double values[] = {1, 7, 4};
  double *data = observations(tree, names, values, 3);
  bt_pic_result *result = NULL;
  assert(bt_pic(tree, data, &result) == CRPX_OK);
  assert(result->count == 2 && result->node[0] == find(tree, "J") && result->node[1] == find(tree, "R"));
  close_to(result->contrast[0], -6 / sqrt(8)); close_to(result->variance[0], 8);
  close_to(result->contrast[1], 0); close_to(result->variance[1], 11);
  bt_pic_free(result);
  data[find(tree, "E")] = NAN;
  assert(bt_pic(tree, data, &result) == CRPX_OK);
  assert(result->count == 1 && result->node[0] == find(tree, "J"));
  close_to(result->contrast[0], -6 / sqrt(8)); close_to(result->variance[0], 8);
  bt_pic_free(result);
  data[find(tree, "C")] = NAN; data[find(tree, "B")] = 4;
  assert(bt_pic(tree, data, &result) == CRPX_OK);
  assert(result->count == 1 && result->node[0] == find(tree, "I"));
  close_to(result->contrast[0], -3 / sqrt(3)); close_to(result->variance[0], 3);
  bt_pic_free(result);
  data[find(tree, "B")] = 1;
  assert(bt_pic(tree, data, &result) == CRPX_OK);
  close_to(result->contrast[0], 0); bt_pic_free(result);
  data[find(tree, "A")] = NAN;
  assert(bt_pic(tree, data, &result) == CRPX_ERROR_INVALID_ARGUMENT && !result);
  data[find(tree, "A")] = INFINITY;
  assert(bt_pic(tree, data, &result) == CRPX_ERROR_INVALID_ARGUMENT && !result);
  free(data); crpx_tree_free(tree);
}

static void missing_regression(void)
{
  crpx_tree *tree = parse("(((A:1,B:2)I:1,C:4)J:2,(D:1,E:1)K:2)R;");
  crpx_tree *observed = parse("(((A:1,B:2)I:1,C:4)J:2)R;");
  const char *names[] = {"A", "B", "C", "D", "E"};
  double response[] = {1, 2, 4, NAN, 10000}, predictor[] = {0, 1, 2, 3, NAN};
  double *y = observations(tree, names, response, 5), *x = observations(tree, names, predictor, 5);
  double *ref_y = observations(observed, names, response, 3), *ref_x = observations(observed, names, predictor, 3);
  bt_pgls_result *result = NULL, *reference = NULL;
  assert(bt_pgls(tree, y, x, 1, &result) == CRPX_OK);
  assert(bt_pgls(observed, ref_y, ref_x, 1, &reference) == CRPX_OK);
  assert(result->tips == 3 && result->parameters == 2);
  close_to(result->sigma2, reference->sigma2); close_to(result->log_likelihood, reference->log_likelihood);
  for (size_t j = 0; j < 2; ++j) {
    close_to(result->coefficients[j], reference->coefficients[j]);
    close_to(result->standard_errors[j], reference->standard_errors[j]);
  }
  bt_pgls_free(result); bt_pgls_free(reference);
  y[find(tree, "B")] = NAN;
  assert(bt_pgls(tree, y, x, 1, &result) == CRPX_ERROR_INVALID_ARGUMENT && !result);
  y[find(tree, "B")] = 2;
  x[find(tree, "D")] = INFINITY;
  assert(bt_pgls(tree, y, x, 1, &result) == CRPX_ERROR_INVALID_ARGUMENT && !result);
  x[find(tree, "D")] = 3;
  for (size_t node = 0; node < crpx_tree_node_count(tree); ++node) y[node] = NAN;
  assert(bt_pgls(tree, y, x, 1, &result) == CRPX_ERROR_INVALID_ARGUMENT && !result);
  free(y); free(x); free(ref_y); free(ref_x); crpx_tree_free(tree); crpx_tree_free(observed);
}

static void independent_columns(void)
{
  crpx_tree *tree = parse("((A:1,B:2)I:1,(C:1,D:3)J:2)R;");
  const char *names[] = {"A", "B", "C", "D"};
  const double first[] = {1, NAN, 4, 7}, second[] = {NAN, 2, NAN, 9};
  size_t nodes = crpx_tree_node_count(tree);
  double *a = observations(tree, names, first, 4), *b = observations(tree, names, second, 4);
  double *matrix = malloc(nodes * 2 * sizeof(*matrix));
  assert(matrix);
  for (size_t node = 0; node < nodes; ++node) {
    matrix[node * 2] = a[node]; matrix[node * 2 + 1] = b[node];
  }
  bt_brownian_batch_result *batch = NULL;
  bt_brownian_result *reference = NULL;
  assert(bt_brownian_fit_many(tree, matrix, 2, &batch) == CRPX_OK);
  assert(batch->traits == 2 && batch->results[0]->tips == 3 && batch->results[1]->tips == 2);
  for (size_t column = 0; column < 2; ++column) {
    assert(bt_brownian_fit(tree, column ? b : a, &reference) == CRPX_OK);
    close_to(batch->results[column]->root_mean, reference->root_mean);
    close_to(batch->results[column]->sigma2, reference->sigma2);
    close_to(batch->results[column]->log_likelihood, reference->log_likelihood);
    for (size_t node = 0; node < nodes; ++node) {
      close_to(batch->results[column]->mean[node], reference->mean[node]);
      close_to(batch->results[column]->variance[node], reference->variance[node]);
    }
    bt_brownian_free(reference);
  }
  bt_brownian_batch_free(batch);
  /* A failed later column does not commit the earlier valid result. */
  matrix[find(tree, "D") * 2 + 1] = NAN;
  batch = (bt_brownian_batch_result *)1;
  assert(bt_brownian_fit_many(tree, matrix, 2, &batch) == CRPX_ERROR_INVALID_ARGUMENT && !batch);
  assert(bt_brownian_fit_many(tree, matrix, 0, &batch) == CRPX_ERROR_INVALID_ARGUMENT && !batch);
  assert(bt_brownian_fit_many(NULL, matrix, 2, &batch) == CRPX_ERROR_INVALID_ARGUMENT && !batch);
  assert(bt_brownian_fit_many(tree, NULL, 2, &batch) == CRPX_ERROR_INVALID_ARGUMENT && !batch);
  assert(bt_brownian_fit_many(tree, matrix, 2, NULL) == CRPX_ERROR_INVALID_ARGUMENT);
  assert(bt_brownian_fit_many(tree, matrix, SIZE_MAX, &batch) == CRPX_ERROR_OVERFLOW && !batch);
  bt_brownian_batch_free(NULL);
  free(matrix); free(a); free(b); crpx_tree_free(tree);
}

static void errors(void)
{
  const char *names[] = {"A", "B"};
  double values[] = {1, 2};
  crpx_tree *tree = parse("(A:0,B:1);");
  double *data = observations(tree, names, values, 2);
  bt_brownian_result *result = (bt_brownian_result *)1;
  assert(bt_brownian_fit(tree, data, &result) == CRPX_ERROR_INVALID_ARGUMENT && !result);
  assert(crpx_tree_set_branch_length(tree, find(tree, "A"), 1.0) == CRPX_OK);
  data[find(tree, "A")] = NAN;
  assert(bt_brownian_fit(tree, data, &result) == CRPX_ERROR_INVALID_ARGUMENT && !result);
  data[find(tree, "A")] = 2.0;
  assert(bt_brownian_fit(tree, data, &result) == CRPX_ERROR_INVALID_ARGUMENT && !result);
  data[find(tree, "A")] = 1.0;
  assert(crpx_tree_clear_branch_length(tree, find(tree, "A")) == CRPX_OK);
  assert(bt_brownian_fit(tree, data, &result) == CRPX_ERROR_NOT_FOUND && !result);
  assert(bt_brownian_fit(NULL, data, &result) == CRPX_ERROR_INVALID_ARGUMENT && !result);
  assert(bt_brownian_fit(tree, data, NULL) == CRPX_ERROR_INVALID_ARGUMENT);
  free(data); crpx_tree_free(tree);
}

int main(void)
{
  brownian_dense_reference(); polytomy_and_unary(); contrasts(); regression(); nested_marginals(); multiple_predictors(); errors();
  missing_brownian(); missing_contrasts(); missing_regression();
  independent_columns();
  puts("Continuous model tests passed");
  return 0;
}
