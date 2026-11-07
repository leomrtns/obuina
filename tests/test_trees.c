/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "obuina.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static crpx_tree *parse(const char *text)
{
  crpx_tree *tree = NULL;
  size_t offset = 0;
  assert(crpx_tree_from_newick(text, &tree, &offset) == CRPX_OK);
  return tree;
}

static crpx_node_id named_tip(const crpx_tree *tree, const char *name)
{
  for (size_t i = 0; i < crpx_tree_node_count(tree); ++i) {
    const char *label = crpx_tree_label(tree, (crpx_node_id)i);
    if (crpx_tree_first_child(tree, (crpx_node_id)i) == CRPX_NODE_NONE && label && !strcmp(name, label))
      return (crpx_node_id)i;
  }
  assert(0);
  return CRPX_NODE_NONE;
}

static double tip_distance(const crpx_tree *tree, const char *a, const char *b)
{
  size_t n = crpx_tree_node_count(tree);
  double *distance = malloc(n * sizeof(*distance));
  assert(distance);
  for (size_t i = 0; i < n; ++i) distance[i] = NAN;
  crpx_node_id node = named_tip(tree, a);
  double total = 0;
  while (node != CRPX_NODE_NONE) {
    distance[node] = total;
    crpx_node_id p = crpx_tree_parent(tree, node);
    if (p != CRPX_NODE_NONE) {
      double edge;
      assert(crpx_tree_branch_length(tree, node, &edge) == CRPX_OK);
      total += edge;
    }
    node = p;
  }
  total = 0;
  node = named_tip(tree, b);
  while (!isfinite(distance[node])) {
    double edge;
    assert(crpx_tree_branch_length(tree, node, &edge) == CRPX_OK);
    total += edge;
    node = crpx_tree_parent(tree, node);
    assert(node != CRPX_NODE_NONE);
  }
  total += distance[node];
  free(distance);
  return total;
}

static double root_depth(const crpx_tree *tree, const char *name)
{
  crpx_node_id node = named_tip(tree, name);
  double total = 0;
  while (crpx_tree_parent(tree, node) != CRPX_NODE_NONE) {
    double edge;
    assert(crpx_tree_branch_length(tree, node, &edge) == CRPX_OK);
    total += edge;
    node = crpx_tree_parent(tree, node);
  }
  return total;
}

static void assert_close(double a, double b)
{
  assert(fabs(a - b) < 1e-10 * fmax(1, fmax(fabs(a), fabs(b))));
}

static void reroot_test(const char *newick, const char *const *names, size_t tips,
                        double expected_height)
{
  crpx_tree *tree = parse(newick), *rooted = NULL;
  assert(bt_tree_midpoint_root(tree, &rooted) == CRPX_OK);
  assert(rooted);
  assert(crpx_tree_leaf_count(rooted) == tips);
  double maximum = 0;
  for (size_t i = 0; i < tips; ++i) {
    for (size_t j = i + 1; j < tips; ++j)
      assert_close(tip_distance(tree, names[i], names[j]), tip_distance(rooted, names[i], names[j]));
    double depth = root_depth(rooted, names[i]);
    if (depth > maximum) maximum = depth;
  }
  assert_close(maximum, expected_height);
  bt_tree_info info;
  assert(bt_tree_inspect(rooted, 1e-10, &info) == CRPX_OK);
  assert(info.tips == tips);
  crpx_tree_free(rooted); crpx_tree_free(tree);
}

static void test_inspect(void)
{
  crpx_tree *tree = parse("((A:1,B:2):3,C:7):99;");
  bt_tree_info info;
  assert(bt_tree_inspect(tree, 1e-10, &info) == CRPX_OK);
  assert(info.nodes == 5 && info.tips == 3 && !info.ultrametric);
  assert_close(info.minimum_tip_depth, 4);
  assert_close(info.maximum_tip_depth, 7);
  assert(bt_tree_inspect(tree, -1, &info) == CRPX_ERROR_INVALID_ARGUMENT);
  assert(bt_tree_inspect(tree, NAN, &info) == CRPX_ERROR_INVALID_ARGUMENT);
  assert(bt_tree_inspect(tree, 0, NULL) == CRPX_ERROR_INVALID_ARGUMENT);
  crpx_tree_free(tree);
  tree = parse("(A:1,A:2);");
  assert(bt_tree_inspect(tree, 1e-10, &info) == CRPX_ERROR_INVALID_ARGUMENT);
  crpx_tree_free(tree);
  tree = parse("(A:1,B);");
  assert(bt_tree_inspect(tree, 1e-10, &info) == CRPX_ERROR_NOT_FOUND);
  crpx_tree_free(tree);
  tree = parse("(A:-1,B:2);");
  assert(bt_tree_inspect(tree, 1e-10, &info) == CRPX_ERROR_INVALID_ARGUMENT);
  crpx_tree_free(tree);
  tree = parse("(:1,B:2);");
  assert(bt_tree_inspect(tree, 1e-10, &info) == CRPX_ERROR_INVALID_ARGUMENT);
  crpx_tree_free(tree);
  tree = parse("(A:0,B:2):-7;");
  assert(bt_tree_inspect(tree, 1e-10, &info) == CRPX_OK); /* Stem is deliberately excluded. */
  assert_close(info.minimum_tip_depth, 0);
  assert_close(info.maximum_tip_depth, 2);
  crpx_tree_free(tree);
}

static void test_midpoint(void)
{
  const char *abc[] = { "A", "B", "C" };
  reroot_test("(A:1,B:2,C:7);", abc, 3, 4.5); /* Midpoint inside edge to C. */
  reroot_test("(A:3,B:3,C:1);", abc, 3, 3); /* Exact original vertex. */
  reroot_test("((A:2,B:1):1,C:4);", abc, 3, 3.5); /* Suppress artificial old root. */
  reroot_test("(((A:2,B:1):1,C:4):9);", abc, 3, 3.5); /* Remove unary stem. */
  reroot_test("(A:0,B:0,C:0);", abc, 3, 0); /* Unidentifiable: clone. */
  reroot_test("((A:0,B:0):2,C:2);", abc, 3, 2);
  const char *ab[] = { "A", "B" };
  reroot_test("(A:1,B:3);", ab, 2, 2);
  reroot_test("(A:2,B:2);", ab, 2, 2);
  const char *abcd[] = { "A", "B", "C", "D" };
  reroot_test("(((A:3,B:2):1):1,(C:1,D:1):1);", abcd, 4, 3.5);
  crpx_tree *tree = parse("(A:1);"), *rooted = (crpx_tree *)(uintptr_t)1;
  assert(bt_tree_midpoint_root(tree, &rooted) == CRPX_ERROR_UNSUPPORTED);
  assert(rooted == NULL);
  crpx_tree_free(tree);
}

static void test_wide_midpoint(void)
{
  /* Construction and midpoint rebuilding must use upstream's tail-aware builder:
   * a wide star exposes the sibling-scan path that costs quadratically. Avoid
   * elapsed-time thresholds; upstream separately instruments constant-time appends. */
  const size_t tips = 4096;
  crpx_tree *tree = NULL, *rooted = NULL;
  crpx_tree_builder *builder = NULL;
  crpx_node_id root, node;
  assert(crpx_tree_new(&tree) == CRPX_OK);
  assert(crpx_tree_builder_new(&builder) == CRPX_OK);
  assert(crpx_tree_reserve(tree, tips + 1) == CRPX_OK);
  assert(crpx_tree_builder_append(builder, tree, CRPX_NODE_NONE, &root) == CRPX_OK);
  assert(crpx_tree_set_label(tree, root, "centre") == CRPX_OK);
  for (size_t i = 0; i < tips; ++i) {
    char name[32];
    snprintf(name, sizeof(name), "tip%zu", i);
    assert(crpx_tree_builder_append(builder, tree, root, &node) == CRPX_OK);
    assert(crpx_tree_set_label(tree, node, name) == CRPX_OK);
    double length = i + 1 == tips ? 9 : 1 + (double)i / (double)tips;
    assert(crpx_tree_set_branch_length(tree, node, length) == CRPX_OK);
  }
  crpx_tree_builder_free(builder); /* Release the borrowed identity before freeing tree. */
  assert(bt_tree_midpoint_root(tree, &rooted) == CRPX_OK);
  assert(crpx_tree_node_count(rooted) == tips + 2);
  assert(crpx_tree_leaf_count(rooted) == tips);
  assert(crpx_tree_validate(rooted) == CRPX_OK);
  const char *names[] = { "tip0", "tip1", "tip2048", "tip4094", "tip4095" };
  for (size_t i = 0; i < sizeof(names) / sizeof(*names); ++i) {
    for (size_t j = i + 1; j < sizeof(names) / sizeof(*names); ++j)
      assert_close(tip_distance(tree, names[i], names[j]), tip_distance(rooted, names[i], names[j]));
  }
  double diameter = 9 + 1 + (double)(tips - 2) / (double)tips;
  assert_close(root_depth(rooted, "tip4095"), diameter / 2);
  assert_close(root_depth(rooted, "tip4094"), diameter / 2);
  size_t centre_count = 0;
  for (size_t i = 0; i < crpx_tree_node_count(rooted); ++i) {
    const char *label = crpx_tree_label(rooted, (crpx_node_id)i);
    if (label && !strcmp(label, "centre")) ++centre_count;
  }
  assert(centre_count == 1);
  assert_close(root_depth(tree, "tip4095"), 9); /* Source remains unchanged. */
  crpx_tree_free(rooted); crpx_tree_free(tree);
}

static void test_extension(void)
{
  crpx_tree *tree = parse("((A:1,B:2):3,C:7);");
  crpx_tree *extended = NULL;
  double added = -1;
  assert(bt_tree_extend_tips(tree, &extended, &added) == CRPX_OK);
  assert_close(added, 5);
  bt_tree_info info;
  assert(bt_tree_inspect(extended, 1e-12, &info) == CRPX_OK && info.ultrametric);
  assert_close(root_depth(extended, "A"), 7);
  assert_close(root_depth(extended, "B"), 7);
  assert_close(root_depth(extended, "C"), 7);
  assert_close(root_depth(tree, "A"), 4); /* Input untouched. */
  crpx_tree_free(extended); crpx_tree_free(tree);
}

static void test_upgma(void)
{
  const char *names[] = { "A", "B", "C" };
  double data[] = { 0, 0, 2, 0, 10, 0 };
  crpx_tree *tree = NULL;
  assert(bt_tree_upgma(names, 3, 2, data, &tree) == CRPX_OK);
  bt_tree_info info;
  assert(bt_tree_inspect(tree, 1e-12, &info) == CRPX_OK && info.ultrametric);
  assert_close(tip_distance(tree, "A", "B"), 2);
  assert_close(tip_distance(tree, "A", "C"), 9);
  assert_close(tip_distance(tree, "B", "C"), 9);
  crpx_tree_free(tree);
  const char *duplicates[] = { "A", "A", "B" };
  assert(bt_tree_upgma(duplicates, 3, 2, data, &tree) == CRPX_ERROR_INVALID_ARGUMENT && !tree);
  data[1] = NAN;
  assert(bt_tree_upgma(names, 3, 2, data, &tree) == CRPX_ERROR_INVALID_ARGUMENT && !tree);
  data[1] = 0;
  assert(bt_tree_upgma(names, SIZE_MAX, 2, data, &tree) == CRPX_ERROR_OVERFLOW && !tree);
  assert(bt_tree_upgma(names, 1, 2, data, &tree) == CRPX_ERROR_INVALID_ARGUMENT && !tree);
}

int main(void)
{
  test_inspect(); test_midpoint(); test_wide_midpoint(); test_extension(); test_upgma();
  puts("tree preparation tests passed");
  return 0;
}
