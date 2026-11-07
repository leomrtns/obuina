/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "obuina.h"
#include <continuous_traits.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define BT_LOG_2PI 1.837877066409345483560659472811

typedef struct {
  size_t nodes, tips;
  crpx_node_id root;
  crpx_tree_workspace *workspace;
  double *length;
} continuous_tree;

static void tree_clear(continuous_tree *index)
{
  crpx_tree_workspace_free(index->workspace);
  free(index->length);
}

/* NaN is an unobserved continuous trait. Infinity remains an invalid value. */
static crpx_status tree_prepare(const crpx_tree *tree, const double *values, continuous_tree *index)
{
  if (!tree || !values) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_tree_validate_policy(tree,
    CRPX_TREE_REQUIRE_COMPLETE_LENGTHS | CRPX_TREE_REQUIRE_POSITIVE_LENGTHS);
  if (status != CRPX_OK) return status;
  index->nodes = crpx_tree_node_count(tree);
  index->tips = crpx_tree_leaf_count(tree);
  index->root = crpx_tree_root(tree);
  if (!index->nodes || index->tips < 2) return CRPX_ERROR_TOPOLOGY;
  if (index->nodes > SIZE_MAX / sizeof(double)) return CRPX_ERROR_OVERFLOW;
  index->length = calloc(index->nodes, sizeof(double));
  if (!index->length) return CRPX_ERROR_ALLOCATION;
  for (size_t node = 0; node < index->nodes; ++node) {
    if (node != index->root) {
      status = crpx_tree_branch_length(tree, (crpx_node_id)node, &index->length[node]);
      if (status != CRPX_OK) return status;
    }
    if (crpx_tree_first_child(tree, (crpx_node_id)node) == CRPX_NODE_NONE && isinf(values[node]))
      return CRPX_ERROR_INVALID_ARGUMENT;
  }
  status = crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_TRAVERSAL | CRPX_TREE_WORKSPACE_DEPTH,
                                 &index->workspace);
  if (status != CRPX_OK) return status;
  return crpx_tree_workspace_refresh(index->workspace, tree);
}

typedef struct {
  const crpx_tree *tree;
  crpx_tree *owned;
  size_t tips;
  crpx_node_id *old_to_new, *new_to_old, *tip_nodes;
  double *values;
} observed_tree;

static void observed_clear(observed_tree *view)
{
  crpx_tree_free(view->owned);
  free(view->old_to_new); free(view->new_to_old); free(view->tip_nodes); free(view->values);
}

/* Missing tips integrate to one in the likelihood. Retain their observed-tip
 * ancestor paths, including the original root for Brownian likelihoods. For PIC
 * contract unary paths and start at the observed MRCA, as contrasts remove the
 * common ancestral state. IDs map results back to the unchanged original tree. */
static crpx_status observed_prepare(const crpx_tree *tree, const double *values,
                                    const continuous_tree *index, bool contract, observed_tree *view)
{
  size_t n = index->nodes;
  if (n > SIZE_MAX / sizeof(crpx_node_id) || n > SIZE_MAX / sizeof(size_t)) return CRPX_ERROR_OVERFLOW;
  view->old_to_new = malloc(n * sizeof(crpx_node_id));
  view->new_to_old = malloc(n * sizeof(crpx_node_id));
  view->tip_nodes = malloc(index->tips * sizeof(crpx_node_id));
  view->values = malloc(index->tips * sizeof(double));
  if (!view->old_to_new || !view->new_to_old || !view->tip_nodes || !view->values)
    return CRPX_ERROR_ALLOCATION;
  for (size_t node = 0; node < n; ++node) {
    view->old_to_new[node] = CRPX_NODE_NONE;
    if (crpx_tree_first_child(tree, (crpx_node_id)node) == CRPX_NODE_NONE && isfinite(values[node]))
      ++view->tips;
  }
  if (view->tips < 2) return CRPX_ERROR_INVALID_ARGUMENT;
  if (view->tips == index->tips) {
    view->tree = tree;
    size_t tip = 0;
    for (size_t node = 0; node < n; ++node) {
      view->old_to_new[node] = view->new_to_old[node] = (crpx_node_id)node;
      if (crpx_tree_first_child(tree, (crpx_node_id)node) == CRPX_NODE_NONE) {
        view->tip_nodes[tip] = (crpx_node_id)node; view->values[tip++] = values[node];
      }
    }
    return CRPX_OK;
  }
  size_t *descendants = calloc(n, sizeof(size_t));
  crpx_node_id *nearest = malloc(n * sizeof(crpx_node_id));
  double *distance = calloc(n, sizeof(double));
  crpx_tree_builder *builder = NULL;
  crpx_status status = CRPX_ERROR_ALLOCATION;
  if (!descendants || !nearest || !distance) goto done;
  const crpx_node_id *postorder = crpx_tree_workspace_postorder(index->workspace);
  const crpx_node_id *preorder = crpx_tree_workspace_preorder(index->workspace);
  for (size_t i = 0; i < n; ++i) {
    crpx_node_id node = postorder[i], child = crpx_tree_first_child(tree, node);
    if (child == CRPX_NODE_NONE) descendants[node] = isfinite(values[node]) ? 1 : 0;
    else for (; child != CRPX_NODE_NONE; child = crpx_tree_next_sibling(tree, child))
      descendants[node] += descendants[child];
  }
  status = crpx_tree_new(&view->owned);
  if (status == CRPX_OK) status = crpx_tree_builder_new(&builder);
  if (status != CRPX_OK) goto done;
  size_t tip = 0;
  for (size_t i = 0; i < n; ++i) {
    crpx_node_id node = preorder[i], parent = crpx_tree_parent(tree, node);
    if (!descendants[node]) continue;
    crpx_node_id target_parent = parent == CRPX_NODE_NONE ? CRPX_NODE_NONE : nearest[parent];
    double edge = parent == CRPX_NODE_NONE ? 0 : distance[parent] + index->length[node];
    if (!isfinite(edge)) { status = CRPX_ERROR_OVERFLOW; goto done; }
    size_t children = 0;
    for (crpx_node_id child = crpx_tree_first_child(tree, node); child != CRPX_NODE_NONE;
         child = crpx_tree_next_sibling(tree, child)) if (descendants[child]) ++children;
    if (contract && children == 1) {
      nearest[node] = target_parent; distance[node] = edge;
      continue;
    }
    crpx_node_id mapped;
    status = crpx_tree_builder_append(builder, view->owned, target_parent, &mapped);
    if (status != CRPX_OK) goto done;
    if (target_parent != CRPX_NODE_NONE) {
      status = crpx_tree_set_branch_length(view->owned, mapped, edge);
      if (status != CRPX_OK) goto done;
    }
    view->old_to_new[node] = mapped; view->new_to_old[mapped] = node;
    nearest[node] = mapped; distance[node] = 0;
    if (!children) { view->tip_nodes[tip] = mapped; view->values[tip++] = values[node]; }
  }
  view->tree = view->owned;
 done:
  crpx_tree_builder_free(builder);
  free(descendants); free(nearest); free(distance);
  return status;
}

void bt_brownian_free(bt_brownian_result *result)
{
  if (!result) return;
  free(result->mean); free(result->variance); free(result);
}

crpx_status bt_brownian_fit(const crpx_tree *tree, const double *values, bt_brownian_result **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  continuous_tree index = {0}; observed_tree view = {0};
  crpx_brownian_result *upstream = NULL;
  bt_brownian_result *result = NULL;
  crpx_status status = tree_prepare(tree, values, &index);
  if (status != CRPX_OK) goto done;
  status = observed_prepare(tree, values, &index, false, &view);
  if (status != CRPX_OK) goto done;
  crpx_brownian_options options;
  crpx_brownian_options_init(&options); options.reconstruct = true;
  status = crpx_brownian_fit(view.tree, view.tip_nodes, view.values, view.tips, &options, &upstream);
  if (status != CRPX_OK) goto done;
  result = calloc(1, sizeof(*result));
  if (!result) { status = CRPX_ERROR_ALLOCATION; goto done; }
  result->nodes = index.nodes; result->tips = view.tips;
  result->root_mean = crpx_brownian_result_root_mean(upstream);
  result->sigma2 = crpx_brownian_result_sigma2(upstream);
  result->log_likelihood = crpx_brownian_result_log_likelihood(upstream);
  result->mean = malloc(index.nodes * sizeof(double));
  result->variance = malloc(index.nodes * sizeof(double));
  if (!result->mean || !result->variance) { status = CRPX_ERROR_ALLOCATION; goto done; }
  const double *means = crpx_brownian_result_means(upstream, NULL);
  const double *variances = crpx_brownian_result_variances(upstream, NULL);
  const crpx_node_id *preorder = crpx_tree_workspace_preorder(index.workspace);
  for (size_t i = 0; i < index.nodes; ++i) {
    crpx_node_id node = preorder[i], mapped = view.old_to_new[node];
    if (mapped != CRPX_NODE_NONE) {
      result->mean[node] = means[mapped]; result->variance[node] = variances[mapped];
    } else {
      crpx_node_id parent = crpx_tree_parent(tree, node);
      result->mean[node] = result->mean[parent];
      result->variance[node] = result->variance[parent] + result->sigma2 * index.length[node];
      if (!isfinite(result->variance[node])) { status = CRPX_ERROR_OVERFLOW; goto done; }
    }
  }
  *out = result; result = NULL;
 done:
  crpx_brownian_result_free(upstream); bt_brownian_free(result);
  observed_clear(&view); tree_clear(&index);
  return status;
}

void bt_pic_free(bt_pic_result *result)
{
  if (!result) return;
  free(result->node); free(result->contrast); free(result->variance); free(result);
}

crpx_status bt_pic(const crpx_tree *tree, const double *values, bt_pic_result **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  continuous_tree index = {0}; observed_tree view = {0};
  crpx_pic_result *upstream = NULL;
  bt_pic_result *result = NULL;
  crpx_status status = tree_prepare(tree, values, &index);
  if (status != CRPX_OK) goto done;
  status = crpx_tree_validate_policy(tree, CRPX_TREE_REQUIRE_BINARY);
  if (status == CRPX_ERROR_UNSUPPORTED) status = CRPX_ERROR_TOPOLOGY;
  if (status != CRPX_OK) goto done;
  status = observed_prepare(tree, values, &index, true, &view);
  if (status != CRPX_OK) goto done;
  status = crpx_pic_calculate(view.tree, view.tip_nodes, view.values, view.tips, &upstream);
  if (status != CRPX_OK) goto done;
  result = calloc(1, sizeof(*result));
  if (!result) { status = CRPX_ERROR_ALLOCATION; goto done; }
  result->count = crpx_pic_result_count(upstream);
  result->node = malloc(result->count * sizeof(crpx_node_id));
  result->contrast = malloc(result->count * sizeof(double));
  result->variance = malloc(result->count * sizeof(double));
  if (!result->node || !result->contrast || !result->variance) { status = CRPX_ERROR_ALLOCATION; goto done; }
  const crpx_node_id *nodes = crpx_pic_result_nodes(upstream);
  const double *contrasts = crpx_pic_result_contrasts(upstream), *variances = crpx_pic_result_variances(upstream);
  for (size_t i = 0; i < result->count; ++i) {
    result->node[i] = view.new_to_old[nodes[i]];
    result->contrast[i] = contrasts[i]; result->variance[i] = variances[i];
  }
  *out = result; result = NULL;
 done:
  crpx_pic_result_free(upstream); bt_pic_free(result);
  observed_clear(&view); tree_clear(&index);
  return status;
}

void bt_pgls_free(bt_pgls_result *result)
{
  if (!result) return;
  free(result->coefficients); free(result->standard_errors); free(result);
}

static crpx_status cholesky(double *matrix, size_t n, double *logdet)
{
  *logdet = 0.0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      double value = matrix[i * n + j];
      for (size_t k = 0; k < j; ++k) value -= matrix[i * n + k] * matrix[j * n + k];
      if (!isfinite(value)) return CRPX_ERROR_OVERFLOW;
      if (i == j) {
        if (!(value > 0.0)) return CRPX_ERROR_INVALID_ARGUMENT;
        matrix[i * n + i] = sqrt(value);
        *logdet += log(value);
      } else matrix[i * n + j] = value / matrix[j * n + j];
    }
  }
  return CRPX_OK;
}

crpx_status bt_pgls(const crpx_tree *tree, const double *response, const double *predictors,
                   size_t columns, bt_pgls_result **out)
{
  continuous_tree index = {0};
  bt_pgls_result *result = NULL;
  double *depth = NULL, *covariance = NULL, *design = NULL, *y = NULL, *r = NULL, *qty = NULL, *work = NULL;
  crpx_node_id *tips = NULL;
  size_t *first_tip = NULL, *after_tip = NULL;
  unsigned char *complete = NULL;
  crpx_status status;
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (columns && !predictors) return CRPX_ERROR_INVALID_ARGUMENT;
  if (columns == SIZE_MAX) return CRPX_ERROR_OVERFLOW;
  status = tree_prepare(tree, response, &index);
  if (status != CRPX_OK) goto done;
  size_t n = 0, p = columns + 1;
  if (columns && index.nodes > SIZE_MAX / columns) { status = CRPX_ERROR_OVERFLOW; goto done; }
  complete = calloc(index.nodes, sizeof(*complete));
  if (!complete) { status = CRPX_ERROR_ALLOCATION; goto done; }
  for (size_t node = 0; node < index.nodes; ++node) {
    if (crpx_tree_first_child(tree, (crpx_node_id)node) != CRPX_NODE_NONE) continue;
    bool observed = isfinite(response[node]);
    for (size_t j = 0; j < columns; ++j) {
      double value = predictors[node * columns + j];
      if (isinf(value)) { status = CRPX_ERROR_INVALID_ARGUMENT; goto done; }
      if (isnan(value)) observed = false;
    }
    if (observed) { complete[node] = 1; ++n; }
  }
  if (n <= p) { status = CRPX_ERROR_INVALID_ARGUMENT; goto done; }
  if (n > SIZE_MAX / n / sizeof(double) || p > SIZE_MAX / p / sizeof(double) ||
      n > SIZE_MAX / p / sizeof(double) || (columns && index.nodes > SIZE_MAX / columns)) {
    status = CRPX_ERROR_OVERFLOW; goto done;
  }
  depth = calloc(index.nodes, sizeof(double));
  first_tip = malloc(index.nodes * sizeof(size_t));
  after_tip = malloc(index.nodes * sizeof(size_t));
  tips = malloc(n * sizeof(*tips));
  covariance = malloc(n * n * sizeof(double));
  design = malloc(n * p * sizeof(double));
  y = malloc(n * sizeof(double));
  r = calloc(p * p, sizeof(double));
  qty = calloc(p, sizeof(double));
  work = malloc(p * sizeof(double));
  result = calloc(1, sizeof(*result));
  if (!depth || !first_tip || !after_tip || !tips || !covariance || !design || !y || !r || !qty || !work || !result) {
    status = CRPX_ERROR_ALLOCATION; goto done;
  }
  result->tips = n; result->parameters = p;
  result->coefficients = calloc(p, sizeof(double));
  result->standard_errors = calloc(p, sizeof(double));
  if (!result->coefficients || !result->standard_errors) { status = CRPX_ERROR_ALLOCATION; goto done; }
  const crpx_node_id *preorder = crpx_tree_workspace_preorder(index.workspace);
  size_t tip_count = 0;
  for (size_t i = 0; i < index.nodes; ++i) {
    crpx_node_id node = preorder[i];
    if (node != index.root) depth[node] = depth[crpx_tree_parent(tree, node)] + index.length[node];
    if (!isfinite(depth[node])) { status = CRPX_ERROR_OVERFLOW; goto done; }
    if (crpx_tree_first_child(tree, node) == CRPX_NODE_NONE) {
      first_tip[node] = tip_count;
      if (complete[node]) tips[tip_count++] = node;
      after_tip[node] = tip_count;
    }
  }
  for (size_t i = 0; i < n; ++i) {
    crpx_node_id a = tips[i];
    y[i] = response[a];
    design[i * p] = 1.0;
    for (size_t j = 1; j < p; ++j) {
      double value = predictors[(size_t)a * columns + j - 1];
      if (!isfinite(value)) { status = CRPX_ERROR_INVALID_ARGUMENT; goto done; }
      design[i * p + j] = value;
    }
    covariance[i * n + i] = depth[a];
  }
  /* Leaf intervals are contiguous in preorder. Fill each off-diagonal pair once
   * at its LCA, including trees with long unary paths, in O(nodes + tips^2). */
  const crpx_node_id *postorder = crpx_tree_workspace_postorder(index.workspace);
  for (size_t position = 0; position < index.nodes; ++position) {
    crpx_node_id node = postorder[position], child = crpx_tree_first_child(tree, node);
    if (child == CRPX_NODE_NONE) continue;
    first_tip[node] = first_tip[child];
    size_t seen = after_tip[child];
    for (child = crpx_tree_next_sibling(tree, child); child != CRPX_NODE_NONE;
         child = crpx_tree_next_sibling(tree, child)) {
      for (size_t i = first_tip[node]; i < seen; ++i)
        for (size_t j = first_tip[child]; j < after_tip[child]; ++j)
          covariance[i * n + j] = covariance[j * n + i] = depth[node];
      seen = after_tip[child];
    }
    after_tip[node] = seen;
  }
  double logdet;
  status = cholesky(covariance, n, &logdet);
  if (status != CRPX_OK) goto done;
  /* Whiten y and X by forward substitution with the Brownian Cholesky factor. */
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < i; ++k) {
      y[i] -= covariance[i * n + k] * y[k];
      for (size_t j = 0; j < p; ++j) design[i * p + j] -= covariance[i * n + k] * design[k * p + j];
    }
    y[i] /= covariance[i * n + i];
    for (size_t j = 0; j < p; ++j) design[i * p + j] /= covariance[i * n + i];
  }
  /* Twice-reorthogonalised QR avoids solving the squared-condition-number
   * normal equations. Full rank is checked relative to each input column. */
  for (size_t j = 0; j < p; ++j) {
    double original_norm = 0.0;
    for (size_t i = 0; i < n; ++i) original_norm = hypot(original_norm, design[i * p + j]);
    for (int pass = 0; pass < 2; ++pass) {
      for (size_t k = 0; k < j; ++k) {
        double product = 0.0;
        for (size_t i = 0; i < n; ++i) product += design[i * p + k] * design[i * p + j];
        r[k * p + j] += product;
        for (size_t i = 0; i < n; ++i) design[i * p + j] -= product * design[i * p + k];
      }
    }
    double norm = 0.0;
    for (size_t i = 0; i < n; ++i) norm = hypot(norm, design[i * p + j]);
    if (!isfinite(norm) || !isfinite(original_norm)) { status = CRPX_ERROR_OVERFLOW; goto done; }
    if (!(norm > DBL_EPSILON * (double)n * original_norm)) {
      status = CRPX_ERROR_INVALID_ARGUMENT; goto done;
    }
    r[j * p + j] = norm;
    for (size_t i = 0; i < n; ++i) {
      design[i * p + j] /= norm;
      qty[j] += design[i * p + j] * y[i];
    }
  }
  for (size_t j = p; j-- > 0;) {
    double value = qty[j];
    for (size_t k = j + 1; k < p; ++k) value -= r[j * p + k] * result->coefficients[k];
    result->coefficients[j] = value / r[j * p + j];
    if (!isfinite(result->coefficients[j])) { status = CRPX_ERROR_OVERFLOW; goto done; }
  }
  double sse = 0.0;
  for (size_t i = 0; i < n; ++i) {
    double residual = y[i];
    for (size_t j = 0; j < p; ++j) residual -= design[i * p + j] * qty[j];
    sse += residual * residual;
  }
  if (!isfinite(sse)) { status = CRPX_ERROR_OVERFLOW; goto done; }
  result->sigma2 = sse / (double)n;
  if (!(result->sigma2 > 0.0)) { status = CRPX_ERROR_INVALID_ARGUMENT; goto done; }
  result->log_likelihood = -0.5 * ((double)n * (BT_LOG_2PI + log(result->sigma2) + 1.0) + logdet);
  for (size_t column = 0; column < p; ++column) {
    for (size_t j = 0; j < p; ++j) work[j] = j == column ? 1.0 : 0.0;
    for (size_t j = p; j-- > 0;) {
      for (size_t k = j + 1; k < p; ++k) work[j] -= r[j * p + k] * work[k];
      work[j] /= r[j * p + j];
      result->standard_errors[j] += work[j] * work[j];
    }
  }
  for (size_t j = 0; j < p; ++j) {
    result->standard_errors[j] = sqrt(result->standard_errors[j] * sse / (double)(n - p));
    if (!isfinite(result->standard_errors[j])) { status = CRPX_ERROR_OVERFLOW; goto done; }
  }
  if (!isfinite(result->log_likelihood)) { status = CRPX_ERROR_OVERFLOW; goto done; }
  *out = result; result = NULL;
  status = CRPX_OK;
done:
  free(depth); free(first_tip); free(after_tip); free(tips);
  free(complete);
  free(covariance); free(design); free(y); free(r); free(qty); free(work);
  bt_pgls_free(result); tree_clear(&index);
  return status;
}
