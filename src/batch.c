/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "obuina.h"
#include <stdint.h>
#include <stdlib.h>

void bt_brownian_batch_free(bt_brownian_batch_result *result)
{
  if (!result) return;
  for (size_t i = 0; i < result->traits; ++i) bt_brownian_free(result->results[i]);
  free(result->results);
  free(result);
}

crpx_status bt_brownian_fit_many(const crpx_tree *tree, const double *values, size_t traits,
                                 bt_brownian_batch_result **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!tree || !values || !traits) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t nodes = crpx_tree_node_count(tree);
  if (!nodes) return CRPX_ERROR_TOPOLOGY;
  if (nodes > SIZE_MAX / traits / sizeof(double) || traits > SIZE_MAX / sizeof(void *))
    return CRPX_ERROR_OVERFLOW;
  bt_brownian_batch_result *result = calloc(1, sizeof(*result));
  double *column = malloc(nodes * sizeof(*column));
  if (!result || !column) { free(result); free(column); return CRPX_ERROR_ALLOCATION; }
  result->results = calloc(traits, sizeof(*result->results));
  if (!result->results) { free(result); free(column); return CRPX_ERROR_ALLOCATION; }
  result->traits = traits;
  crpx_status status = CRPX_OK;
  for (size_t j = 0; j < traits; ++j) {
    for (size_t i = 0; i < nodes; ++i)
      column[i] = crpx_tree_first_child(tree, (crpx_node_id)i) == CRPX_NODE_NONE ? values[i * traits + j] : 0;
    status = bt_brownian_fit(tree, column, &result->results[j]);
    if (status != CRPX_OK) break;
  }
  free(column);
  if (status != CRPX_OK) bt_brownian_batch_free(result);
  else *out = result;
  return status;
}
