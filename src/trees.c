/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "obuina.h"
#include <distance_matrix.h>
#include <upgma.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int compare_names(const void *a, const void *b)
{
  return strcmp(*(const char *const *)a, *(const char *const *)b);
}

static crpx_status unique_names(const char *const *names, size_t count)
{
  if (!names || !count) return CRPX_ERROR_INVALID_ARGUMENT;
  if (count > SIZE_MAX / sizeof(char *)) return CRPX_ERROR_OVERFLOW;
  const char **sorted = malloc(count * sizeof(*sorted));
  if (!sorted) return CRPX_ERROR_ALLOCATION;
  crpx_status status = CRPX_OK;
  for (size_t i = 0; i < count; ++i) {
    if (!names[i] || !names[i][0]) { status = CRPX_ERROR_INVALID_ARGUMENT; goto done; }
    sorted[i] = names[i];
  }
  qsort(sorted, count, sizeof(*sorted), compare_names);
  for (size_t i = 1; i < count; ++i)
    if (!strcmp(sorted[i - 1], sorted[i])) { status = CRPX_ERROR_INVALID_ARGUMENT; break; }
 done:
  free(sorted);
  return status;
}

/* Apply this application's requirements through curupixa's shared validator.
 * Root stems remain excluded from evolutionary depth and sign checks. */
static crpx_status inspect(const crpx_tree *tree, double tolerance, bt_tree_info *out, double *depth)
{
  if (!tree || !out || !isfinite(tolerance) || tolerance < 0) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = crpx_tree_validate_policy(tree, CRPX_TREE_REQUIRE_UNIQUE_TIPS |
    CRPX_TREE_REQUIRE_COMPLETE_LENGTHS | CRPX_TREE_REQUIRE_NONNEGATIVE_LENGTHS);
  if (status != CRPX_OK) return status;
  size_t n = crpx_tree_node_count(tree), tips = crpx_tree_leaf_count(tree);
  if (!n || !tips) return CRPX_ERROR_TOPOLOGY;
  if (n > SIZE_MAX / sizeof(double)) return CRPX_ERROR_OVERFLOW;
  double *owned_depth = depth ? NULL : calloc(n, sizeof(*owned_depth));
  crpx_tree_workspace *workspace = NULL;
  if (!depth && !owned_depth) { status = CRPX_ERROR_ALLOCATION; goto done; }
  if (!depth) depth = owned_depth;
  status = crpx_tree_workspace_new(CRPX_TREE_WORKSPACE_TRAVERSAL, &workspace);
  if (status != CRPX_OK) goto done;
  status = crpx_tree_workspace_refresh(workspace, tree);
  if (status != CRPX_OK) goto done;
  const crpx_node_id *order = crpx_tree_workspace_preorder(workspace);
  crpx_node_id root = crpx_tree_root(tree);
  depth[root] = 0;
  bt_tree_info info = { n, tips, INFINITY, 0, false };
  for (size_t i = 0; i < n; ++i) {
    crpx_node_id node = order[i];
    if (node != root) {
      double edge;
      status = crpx_tree_branch_length(tree, node, &edge);
      if (status != CRPX_OK) goto done;
      depth[node] = depth[crpx_tree_parent(tree, node)] + edge;
      if (!isfinite(depth[node])) { status = CRPX_ERROR_OVERFLOW; goto done; }
    }
    if (crpx_tree_first_child(tree, node) == CRPX_NODE_NONE) {
      if (depth[node] < info.minimum_tip_depth) info.minimum_tip_depth = depth[node];
      if (depth[node] > info.maximum_tip_depth) info.maximum_tip_depth = depth[node];
    }
  }
  info.ultrametric = info.maximum_tip_depth - info.minimum_tip_depth <=
    tolerance * fmax(1.0, info.maximum_tip_depth);
  *out = info;
 done:
  crpx_tree_workspace_free(workspace);
  free(owned_depth);
  return status;
}

crpx_status bt_tree_inspect(const crpx_tree *tree, double tolerance, bt_tree_info *out)
{
  return inspect(tree, tolerance, out, NULL);
}

crpx_status bt_tree_extend_tips(const crpx_tree *tree, crpx_tree **out, double *total_added)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  size_t n = crpx_tree_node_count(tree);
  if (!n || n > SIZE_MAX / sizeof(double)) return CRPX_ERROR_INVALID_ARGUMENT;
  double *depth = malloc(n * sizeof(*depth));
  if (!depth) return CRPX_ERROR_ALLOCATION;
  bt_tree_info info;
  crpx_status status = inspect(tree, 1e-10, &info, depth);
  crpx_tree *copy = NULL;
  double added = 0;
  if (status != CRPX_OK) goto done;
  status = crpx_tree_clone(tree, &copy);
  if (status != CRPX_OK) goto done;
  for (size_t i = 0; i < n; ++i) {
    crpx_node_id node = (crpx_node_id)i;
    if (crpx_tree_first_child(tree, node) != CRPX_NODE_NONE || node == crpx_tree_root(tree)) continue;
    double edge, increment = info.maximum_tip_depth - depth[node];
    status = crpx_tree_branch_length(tree, node, &edge);
    if (status != CRPX_OK) goto done;
    added += increment;
    if (!isfinite(added) || !isfinite(edge + increment)) { status = CRPX_ERROR_OVERFLOW; goto done; }
    status = crpx_tree_set_branch_length(copy, node, edge + increment);
    if (status != CRPX_OK) goto done;
  }
  *out = copy;
  copy = NULL;
  if (total_added) *total_added = added;
 done:
  crpx_tree_free(copy);
  free(depth);
  return status;
}

crpx_status bt_tree_upgma(const char *const *names, size_t taxa, size_t features,
                         const double *metadata, crpx_tree **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!metadata || taxa < 2 || !features) return CRPX_ERROR_INVALID_ARGUMENT;
  if (taxa > SIZE_MAX / features || taxa * features > SIZE_MAX / sizeof(double)) return CRPX_ERROR_OVERFLOW;
  crpx_status status = unique_names(names, taxa);
  if (status != CRPX_OK) return status;
  for (size_t i = 0; i < taxa * features; ++i)
    if (!isfinite(metadata[i])) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_distance_matrix_t matrix = NULL;
  status = crpx_distance_matrix_new(taxa, &matrix);
  if (status != CRPX_OK) return status;
  for (size_t i = 1; i < taxa; ++i) {
    for (size_t j = 0; j < i; ++j) {
      double distance = 0;
      for (size_t k = 0; k < features; ++k)
        distance = hypot(distance, metadata[i * features + k] - metadata[j * features + k]);
      if (!isfinite(distance)) { status = CRPX_ERROR_OVERFLOW; goto done; }
      status = crpx_distance_matrix_set(matrix, i, j, distance);
      if (status != CRPX_OK) goto done;
    }
  }
  status = crpx_tree_upgma(matrix, names, out);
 done:
  crpx_distance_matrix_free(matrix);
  return status;
}

typedef struct {
  size_t n;
  size_t *head, *next, *reverse;
  crpx_node_id *to;
  double *length;
  unsigned char *active, *taxon;
} graph;

static void graph_free(graph *g)
{
  free(g->head); free(g->next); free(g->reverse); free(g->to); free(g->length); free(g->active); free(g->taxon);
}

static crpx_status graph_new(const crpx_tree *tree, graph *g)
{
  size_t n = crpx_tree_node_count(tree);
  if (n > SIZE_MAX / (2 * sizeof(double)) || n > SIZE_MAX / (2 * sizeof(size_t)))
    return CRPX_ERROR_OVERFLOW;
  g->n = n;
  g->head = malloc(n * sizeof(*g->head));
  g->next = malloc(2 * n * sizeof(*g->next));
  g->reverse = malloc(2 * n * sizeof(*g->reverse));
  g->to = malloc(2 * n * sizeof(*g->to));
  g->length = malloc(2 * n * sizeof(*g->length));
  g->active = malloc(n);
  g->taxon = malloc(n);
  if (!g->head || !g->next || !g->reverse || !g->to || !g->length || !g->active || !g->taxon)
    return CRPX_ERROR_ALLOCATION;
  for (size_t i = 0; i < n; ++i) {
    g->head[i] = SIZE_MAX;
    g->active[i] = 1;
    g->taxon[i] = crpx_tree_first_child(tree, (crpx_node_id)i) == CRPX_NODE_NONE;
  }
  for (size_t e = 0; e < 2 * n; ++e) g->reverse[e] = e ^ 1;
  size_t edge = 0;
  crpx_node_id root = crpx_tree_root(tree);
  for (size_t i = 0; i < n; ++i) {
    if (i == root) continue;
    crpx_node_id parent = crpx_tree_parent(tree, (crpx_node_id)i);
    double length;
    crpx_status status = crpx_tree_branch_length(tree, (crpx_node_id)i, &length);
    if (status != CRPX_OK) return status;
    g->to[edge] = parent; g->length[edge] = length; g->next[edge] = g->head[i]; g->head[i] = edge++;
    g->to[edge] = (crpx_node_id)i; g->length[edge] = length;
    g->next[edge] = g->head[parent]; g->head[parent] = edge++;
  }
  /* A chain above the first branching vertex has no taxon meaning. Removing it
   * prevents the old root from becoming an extra leaf after reorientation. */
  crpx_node_id node = root;
  while (!g->taxon[node]) {
    size_t degree = 0;
    crpx_node_id next = CRPX_NODE_NONE;
    for (size_t e = g->head[node]; e != SIZE_MAX; e = g->next[e])
      if (g->active[g->to[e]]) { ++degree; next = g->to[e]; }
    if (degree != 1) break;
    g->active[node] = 0;
    node = next;
  }
  /* Suppress unlabelled degree-two vertices, including an artificial binary
   * old root, by summing its incident lengths. Preserve labelled vertices so
   * annotations remain attached to the same graph vertex. */
  for (size_t i = 0; i < n; ++i) {
    const char *label = crpx_tree_label(tree, (crpx_node_id)i);
    if (!g->active[i] || g->taxon[i] || (label && label[0])) continue;
    size_t first = SIZE_MAX, second = SIZE_MAX, degree = 0;
    for (size_t e = g->head[i]; e != SIZE_MAX; e = g->next[e]) {
      if (!g->active[g->to[e]]) continue;
      if (degree == 0) first = e;
      if (degree == 1) second = e;
      ++degree;
    }
    if (degree != 2) continue;
    double length = g->length[first] + g->length[second];
    if (!isfinite(length)) return CRPX_ERROR_OVERFLOW;
    crpx_node_id u = g->to[first], v = g->to[second];
    size_t reverse_first = g->reverse[first], reverse_second = g->reverse[second];
    g->to[reverse_first] = v; g->length[reverse_first] = length;
    g->to[reverse_second] = u; g->length[reverse_second] = length;
    g->reverse[reverse_first] = reverse_second;
    g->reverse[reverse_second] = reverse_first;
    g->active[i] = 0;
  }
  return CRPX_OK;
}

static crpx_status sweep(const graph *g, crpx_node_id start, crpx_node_id *parent,
                         double *distance, crpx_node_id *queue, crpx_node_id *farthest)
{
  size_t used = 1;
  queue[0] = start;
  parent[start] = CRPX_NODE_NONE;
  distance[start] = 0;
  *farthest = start;
  for (size_t at = 0; at < used; ++at) {
    crpx_node_id node = queue[at];
    if (g->taxon[node] && (distance[node] > distance[*farthest] ||
        (distance[node] == distance[*farthest] && node < *farthest))) *farthest = node;
    for (size_t edge = g->head[node]; edge != SIZE_MAX; edge = g->next[edge]) {
      crpx_node_id child = g->to[edge];
      if (!g->active[child] || child == parent[node]) continue;
      parent[child] = node;
      distance[child] = distance[node] + g->length[edge];
      if (!isfinite(distance[child])) return CRPX_ERROR_OVERFLOW;
      queue[used++] = child;
    }
  }
  return CRPX_OK;
}

static crpx_status copy_label(const crpx_tree *source, crpx_node_id old, crpx_tree *target, crpx_node_id node)
{
  const char *label = crpx_tree_label(source, old);
  return label ? crpx_tree_set_label(target, node, label) : CRPX_OK;
}

/* Build both sides of a split edge, or reorient from an existing vertex.
 * Keep the established midpoint contraction/label policy: upstream edge rooting
 * deliberately preserves every original vertex and does not implement it.
 * The optional upstream builder avoids sibling scans at broad polytomies. */
static crpx_status orient(const crpx_tree *tree, const graph *g, crpx_node_id a, crpx_node_id b,
                          double da, double db, crpx_tree **out)
{
  crpx_tree *copy = NULL;
  crpx_tree_builder *builder = NULL;
  crpx_node_id *queue = malloc(g->n * sizeof(*queue));
  crpx_node_id *parent = malloc(g->n * sizeof(*parent));
  crpx_node_id *map = malloc(g->n * sizeof(*map));
  crpx_status status = CRPX_ERROR_ALLOCATION;
  if (!queue || !parent || !map) goto done;
  status = crpx_tree_new(&copy);
  if (status != CRPX_OK) goto done;
  status = crpx_tree_builder_new(&builder);
  if (status != CRPX_OK) goto done;
  crpx_node_id root;
  status = crpx_tree_builder_append(builder, copy, CRPX_NODE_NONE, &root);
  if (status != CRPX_OK) goto done;
  size_t used = 0;
  if (b == CRPX_NODE_NONE) {
    queue[used++] = a; parent[a] = CRPX_NODE_NONE; map[a] = root;
    status = copy_label(tree, a, copy, root);
    if (status != CRPX_OK) goto done;
  } else {
    queue[used++] = a; queue[used++] = b; parent[a] = b; parent[b] = a;
    status = crpx_tree_builder_append(builder, copy, root, &map[a]);
    if (status != CRPX_OK) goto done;
    status = crpx_tree_set_branch_length(copy, map[a], da);
    if (status != CRPX_OK) goto done;
    status = copy_label(tree, a, copy, map[a]);
    if (status != CRPX_OK) goto done;
    status = crpx_tree_builder_append(builder, copy, root, &map[b]);
    if (status != CRPX_OK) goto done;
    status = crpx_tree_set_branch_length(copy, map[b], db);
    if (status != CRPX_OK) goto done;
    status = copy_label(tree, b, copy, map[b]);
    if (status != CRPX_OK) goto done;
  }
  for (size_t at = 0; at < used; ++at) {
    crpx_node_id node = queue[at];
    for (size_t edge = g->head[node]; edge != SIZE_MAX; edge = g->next[edge]) {
      crpx_node_id child = g->to[edge];
      if (!g->active[child] || child == parent[node]) continue;
      queue[used++] = child; parent[child] = node;
      status = crpx_tree_builder_append(builder, copy, map[node], &map[child]);
      if (status != CRPX_OK) goto done;
      status = crpx_tree_set_branch_length(copy, map[child], g->length[edge]);
      if (status != CRPX_OK) goto done;
      status = copy_label(tree, child, copy, map[child]);
      if (status != CRPX_OK) goto done;
    }
  }
  status = crpx_tree_validate(copy);
  if (status == CRPX_OK) { *out = copy; copy = NULL; }
 done:
  crpx_tree_builder_free(builder);
  crpx_tree_free(copy);
  free(queue); free(parent); free(map);
  return status;
}

crpx_status bt_tree_midpoint_root(const crpx_tree *tree, crpx_tree **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  bt_tree_info info;
  crpx_status status = bt_tree_inspect(tree, 1e-10, &info);
  if (status != CRPX_OK) return status;
  if (info.tips < 2) return CRPX_ERROR_UNSUPPORTED;
  graph g = { 0 };
  crpx_node_id *parent = NULL, *queue = NULL;
  double *distance = NULL;
  status = graph_new(tree, &g);
  if (status != CRPX_OK) goto done;
  parent = malloc(g.n * sizeof(*parent));
  queue = malloc(g.n * sizeof(*queue));
  distance = malloc(g.n * sizeof(*distance));
  if (!parent || !queue || !distance) { status = CRPX_ERROR_ALLOCATION; goto done; }
  crpx_node_id start = 0, a, b;
  while (!g.taxon[start]) ++start;
  status = sweep(&g, start, parent, distance, queue, &a);
  if (status != CRPX_OK) goto done;
  status = sweep(&g, a, parent, distance, queue, &b);
  if (status != CRPX_OK) goto done;
  double diameter = distance[b], midpoint = diameter / 2;
  if (diameter == 0) { status = crpx_tree_clone(tree, out); goto done; }
  crpx_node_id node = b;
  while (parent[node] != CRPX_NODE_NONE && distance[parent[node]] > midpoint) node = parent[node];
  double tolerance = 16 * DBL_EPSILON * diameter;
  if (fabs(distance[node] - midpoint) <= tolerance && !g.taxon[node]) {
    status = orient(tree, &g, node, CRPX_NODE_NONE, 0, 0, out);
  } else {
    crpx_node_id p = parent[node];
    if (p == CRPX_NODE_NONE) { status = CRPX_ERROR_TOPOLOGY; goto done; }
    if (fabs(distance[p] - midpoint) <= tolerance && !g.taxon[p])
      status = orient(tree, &g, p, CRPX_NODE_NONE, 0, 0, out);
    else
      status = orient(tree, &g, node, p, distance[node] - midpoint, midpoint - distance[p], out);
  }
 done:
  free(parent); free(queue); free(distance);
  graph_free(&g);
  return status;
}
