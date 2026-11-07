/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "obuina.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void)
{
  crpx_tree *tree = NULL;
  bt_brownian_result *fit = NULL;
  crpx_status status = crpx_tree_from_newick("((A:1,B:2):1,C:2);", &tree, NULL);
  if (status != CRPX_OK) return 1;
  size_t count = crpx_tree_node_count(tree);
  double *values = malloc(count * sizeof(*values));
  if (!values) { crpx_tree_free(tree); return 1; }
  for (size_t i = 0; i < count; ++i) {
    const char *label = crpx_tree_label(tree, (crpx_node_id)i);
    values[i] = label && !strcmp(label, "A") ? 1 : label && !strcmp(label, "B") ? 3 : 4;
  }
  status = bt_brownian_fit(tree, values, &fit);
  if (status == CRPX_OK)
    printf("root mean = %.6f, diffusion variance = %.6f, log likelihood = %.6f\n",
      fit->root_mean, fit->sigma2, fit->log_likelihood);
  else fprintf(stderr, "%s\n", crpx_status_string(status));
  bt_brownian_free(fit);
  free(values);
  crpx_tree_free(tree);
  return status == CRPX_OK ? 0 : 1;
}
