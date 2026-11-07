/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "obuina.h"
#include <substitution_model.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define OK(call) assert((call) == CRPX_OK)
static void close_value(double a, double b)
{ assert(fabs(a-b) < 2e-9*(1+fabs(b))); }

/* Independently enumerate every joint state assignment, including the observed
 * tips. This checks marginals as well as likelihood and both optimal sets. */
static void exhaustive(const char *newick, size_t k, bool missing)
{
  crpx_tree *tree;
  OK(crpx_tree_from_newick(newick,&tree,NULL));
  size_t n = crpx_tree_node_count(tree), count = 1;
  size_t *observed = malloc(n*sizeof(size_t)), *assignment = malloc(n*sizeof(size_t));
  double *marginal = calloc(n*k,sizeof(double)), *matrix = calloc(n*k*k,sizeof(double));
  unsigned char *optimal = calloc(n*k,sizeof(unsigned char));
  assert(observed && assignment && marginal && matrix && optimal);
  size_t tip = 0;
  for (size_t i = 0; i < n; ++i) {
    count *= k; observed[i] = SIZE_MAX;
    if (crpx_tree_first_child(tree,(crpx_node_id)i) == CRPX_NODE_NONE) {
      observed[i] = missing && tip == 1 ? SIZE_MAX : (2*tip+tip/2)%k; ++tip;
    }
  }
  const double rate = 0.73;
  crpx_substitution_model *model;
  OK(crpx_substitution_model_poisson(k,&model));
  for (size_t i = 0; i < n; ++i) if ((crpx_node_id)i != crpx_tree_root(tree)) {
    double length;
    OK(crpx_tree_branch_length(tree,(crpx_node_id)i,&length));
    OK(crpx_substitution_model_transition(model,rate*length,CRPX_BRANCH_FIXED,matrix+i*k*k,k*k));
  }
  long double total = 0;
  for (size_t code = 0; code < count; ++code) {
    size_t value = code;
    bool allowed = true;
    for (size_t i = 0; i < n; ++i) {
      assignment[i] = value%k; value /= k;
      if (observed[i] != SIZE_MAX && assignment[i] != observed[i]) allowed = false;
    }
    if (!allowed) continue;
    long double probability = 1.0L/k;
    for (size_t i = 0; i < n; ++i) {
      crpx_node_id parent = crpx_tree_parent(tree,(crpx_node_id)i);
      if (parent != CRPX_NODE_NONE) probability *= matrix[(i*k+assignment[parent])*k+assignment[i]];
    }
    total += probability;
    for (size_t i = 0; i < n; ++i) marginal[i*k+assignment[i]] += (double)probability;
  }
  bt_discrete_result *d;
  OK(bt_discrete_asr(tree,observed,k,rate,false,&d));
  close_value(d->log_likelihood,(double)logl(total)); close_value(d->rate,rate);
  for (size_t i = 0; i < n*k; ++i) close_value(d->probabilities[i],marginal[i]/(double)total);
  bt_discrete_free(d);
  for (unsigned ordered = 0; ordered < 2; ++ordered) {
    double minimum = INFINITY;
    memset(optimal,0,n*k);
    for (size_t code = 0; code < count; ++code) {
      size_t value = code;
      bool allowed = true;
      for (size_t i = 0; i < n; ++i) {
        assignment[i] = value%k; value /= k;
        if (observed[i] != SIZE_MAX && assignment[i] != observed[i]) allowed = false;
      }
      if (!allowed) continue;
      double score = 0;
      for (size_t i = 0; i < n; ++i) {
        crpx_node_id parent = crpx_tree_parent(tree,(crpx_node_id)i);
        if (parent != CRPX_NODE_NONE) score += ordered ? fabs((double)assignment[parent]-assignment[i]) :
          (assignment[parent] != assignment[i] ? 1 : 0);
      }
      if (score < minimum) { minimum = score; memset(optimal,0,n*k); }
      if (score == minimum) for (size_t i = 0; i < n; ++i) optimal[i*k+assignment[i]] = 1;
    }
    bt_parsimony_result *p;
    OK(bt_parsimony_asr(tree,observed,k,ordered != 0,&p));
    close_value(p->score,minimum);
    assert(memcmp(p->optimal,optimal,n*k) == 0);
    double reconstructed = 0;
    for (size_t i = 0; i < n; ++i) {
      assert(p->assignment[i] < k && optimal[i*k+p->assignment[i]]);
      if (observed[i] != SIZE_MAX) assert(p->assignment[i] == observed[i]);
      crpx_node_id parent = crpx_tree_parent(tree,(crpx_node_id)i);
      if (parent != CRPX_NODE_NONE)
        reconstructed += ordered ? fabs((double)p->assignment[parent]-p->assignment[i]) :
          (p->assignment[parent] != p->assignment[i] ? 1 : 0);
    }
    close_value(reconstructed,minimum);
    bt_parsimony_free(p);
  }
  crpx_substitution_model_free(model); crpx_tree_free(tree);
  free(observed); free(assignment); free(marginal); free(matrix); free(optimal);
}

static void special_cases(void)
{
  crpx_tree *tree;
  bt_discrete_result *d;
  bt_parsimony_result *p;
  OK(crpx_tree_from_newick("(A:0,B:0);",&tree,NULL));
  size_t observations[] = {SIZE_MAX,0,1};
  assert(bt_discrete_asr(tree,observations,2,1,false,&d) == CRPX_ERROR_INVALID_ARGUMENT && !d);
  OK(bt_parsimony_asr(tree,observations,2,false,&p)); close_value(p->score,1); bt_parsimony_free(p);
  observations[2] = 0;
  OK(bt_discrete_asr(tree,observations,2,1,false,&d));
  close_value(d->log_likelihood,-log(2.0)); close_value(d->probabilities[0],1); bt_discrete_free(d);
  assert(bt_discrete_asr(tree,observations,2,1,true,&d) == CRPX_ERROR_INVALID_ARGUMENT && !d);
  assert(bt_discrete_asr(tree,observations,1,1,false,&d) == CRPX_ERROR_INVALID_ARGUMENT && !d);
  assert(bt_discrete_asr(tree,observations,2,0,false,&d) == CRPX_ERROR_INVALID_ARGUMENT && !d);
  observations[2] = 2;
  assert(bt_discrete_asr(tree,observations,2,1,false,&d) == CRPX_ERROR_OUT_OF_RANGE && !d);
  assert(bt_parsimony_asr(tree,observations,2,false,&p) == CRPX_ERROR_OUT_OF_RANGE && !p);
  observations[1] = observations[2] = 0;
  OK(bt_parsimony_asr(tree,observations,1,false,&p)); close_value(p->score,0); bt_parsimony_free(p);
  crpx_tree_free(tree);
  OK(crpx_tree_from_newick("((A:0.1,B:0.2):0.3,C:0.4);",&tree,NULL));
  size_t states[5] = {SIZE_MAX,SIZE_MAX,0,0,1};
  OK(bt_discrete_asr(tree,states,2,1,true,&d));
  double fitted_ll = d->log_likelihood, fitted_rate = d->rate;
  assert(fitted_rate >= 1e-8/0.25 && fitted_rate <= 100/0.25*(1+1e-12));
  bt_discrete_free(d);
  for (size_t i = 0; i <= 100; ++i) {
    double rate = exp(log(1e-8/0.25)+log(1e10)*(double)i/100);
    OK(bt_discrete_asr(tree,states,2,rate,false,&d));
    assert(fitted_ll+1e-8 >= d->log_likelihood); bt_discrete_free(d);
  }
  for (size_t i = 0; i < 5; ++i) states[i] = SIZE_MAX;
  OK(bt_discrete_asr(tree,states,3,1,false,&d)); close_value(d->log_likelihood,0);
  for (size_t i = 0; i < 15; ++i) close_value(d->probabilities[i],1.0/3);
  bt_discrete_free(d); crpx_tree_free(tree);
  OK(crpx_tree_from_newick("(A:-1,B:0.1);",&tree,NULL));
  assert(bt_discrete_asr(tree,observations,2,1,false,&d) == CRPX_ERROR_INVALID_ARGUMENT && !d);
  OK(bt_parsimony_asr(tree,observations,2,false,&p)); bt_parsimony_free(p); crpx_tree_free(tree);
  OK(crpx_tree_from_newick("(A,B);",&tree,NULL));
  assert(bt_discrete_asr(tree,observations,2,1,false,&d) == CRPX_ERROR_NOT_FOUND && !d);
  crpx_tree_free(tree);
}

static void underflow_case(void)
{
  const size_t taxa = 600;
  crpx_tree *tree;
  crpx_node_id root, tip;
  OK(crpx_tree_new(&tree)); OK(crpx_tree_add_node(tree,CRPX_NODE_NONE,&root));
  size_t *states = malloc((taxa+1)*sizeof(size_t)); assert(states); states[0] = SIZE_MAX;
  for (size_t i = 0; i < taxa; ++i) {
    char label[32]; snprintf(label,sizeof(label),"taxon%zu",i);
    OK(crpx_tree_add_node(tree,root,&tip)); OK(crpx_tree_set_label(tree,tip,label));
    OK(crpx_tree_set_branch_length(tree,tip,1e-200)); states[tip] = i%2;
  }
  bt_discrete_result *d;
  OK(bt_discrete_asr(tree,states,2,1,false,&d));
  close_value(d->log_likelihood,(double)(taxa/2)*log(1e-200));
  close_value(d->probabilities[root*2],0.5);
  for (size_t i = 0; i < d->nodes; ++i)
    close_value(d->probabilities[i*2]+d->probabilities[i*2+1],1);
  bt_discrete_free(d); free(states); crpx_tree_free(tree);
}

static void uncertain_observations(void)
{
  crpx_tree *tree;
  /* Unnamed tips exercise the node-ID adaptation; internal NaNs are ignored. */
  OK(crpx_tree_from_newick("((:0.2,:0.7):0.3,:1.1);", &tree, NULL));
  const size_t n = 5, k = 3;
  const double rate = 0.43;
  double weights[] = { NAN, NAN, NAN, NAN, NAN, NAN,
    0.8, 0.15, 0.05, 1, 1, 1, 0, 2, 1 };
  double expected[15] = {0};
  long double total = 0;
  for (size_t code = 0; code < 243; ++code) {
    size_t assignment[5], value = code;
    for (size_t i = 0; i < n; ++i) { assignment[i] = value % k; value /= k; }
    long double probability = 1.0L / k;
    for (size_t i = 1; i < n; ++i) {
      double length;
      OK(crpx_tree_branch_length(tree, (crpx_node_id)i, &length));
      size_t parent = crpx_tree_parent(tree, (crpx_node_id)i);
      double decay = exp(-(double)k * rate * length / (double)(k - 1));
      probability *= 1.0 / k + ((assignment[parent] == assignment[i] ? 1.0 : 0.0) - 1.0 / k) * decay;
      if (crpx_tree_first_child(tree, (crpx_node_id)i) == CRPX_NODE_NONE)
        probability *= weights[i * k + assignment[i]];
    }
    total += probability;
    for (size_t i = 0; i < n; ++i) expected[i * k + assignment[i]] += (double)probability;
  }
  bt_discrete_result *result, *scaled, *fitted;
  OK(bt_discrete_asr_likelihoods(tree, weights, k, rate, false, &result));
  close_value(result->log_likelihood, (double)logl(total));
  for (size_t i = 0; i < n * k; ++i) close_value(result->probabilities[i], expected[i] / (double)total);
  for (size_t s = 0; s < k; ++s) weights[2 * k + s] *= 7;
  OK(bt_discrete_asr_likelihoods(tree, weights, k, rate, false, &scaled));
  close_value(scaled->log_likelihood, result->log_likelihood + log(7.0));
  for (size_t i = 0; i < n * k; ++i) close_value(scaled->probabilities[i], result->probabilities[i]);
  bt_discrete_free(result); bt_discrete_free(scaled);
  OK(bt_discrete_asr_likelihoods(tree, weights, k, rate, true, &fitted));
  for (size_t step = 0; step < 30; ++step) {
    double candidate = exp(log(1e-8 / 0.575) + log(1e10) * (double)step / 29);
    OK(bt_discrete_asr_likelihoods(tree, weights, k, candidate, false, &result));
    assert(fitted->log_likelihood + 1e-8 >= result->log_likelihood);
    bt_discrete_free(result);
  }
  bt_discrete_free(fitted);
  for (size_t i = 0; i < n; ++i) assert(crpx_tree_label(tree, (crpx_node_id)i) == NULL);

  const double invalid[] = {-1, NAN, INFINITY};
  for (size_t i = 0; i < sizeof(invalid) / sizeof(*invalid); ++i) {
    double old = weights[6]; weights[6] = invalid[i];
    result = (bt_discrete_result *)1;
    assert(bt_discrete_asr_likelihoods(tree, weights, k, rate, false, &result) == CRPX_ERROR_INVALID_ARGUMENT);
    assert(result == NULL); weights[6] = old;
  }
  weights[6] = weights[7] = weights[8] = 0;
  assert(bt_discrete_asr_likelihoods(tree, weights, k, rate, false, &result) == CRPX_ERROR_INVALID_ARGUMENT);
  assert(result == NULL);
  /* Exact, unknown and likelihood-vector input agree, including missing tips. */
  size_t observed[] = {SIZE_MAX, SIZE_MAX, 1, SIZE_MAX, 0};
  for (size_t i = 2; i < n; ++i)
    for (size_t s = 0; s < k; ++s) weights[i * k + s] = observed[i] == SIZE_MAX || observed[i] == s ? 1 : 0;
  OK(bt_discrete_asr(tree, observed, k, rate, false, &result));
  OK(bt_discrete_asr_likelihoods(tree, weights, k, rate, false, &scaled));
  close_value(result->log_likelihood, scaled->log_likelihood);
  for (size_t i = 0; i < n * k; ++i) close_value(result->probabilities[i], scaled->probabilities[i]);
  bt_discrete_free(result); bt_discrete_free(scaled);
  crpx_tree_free(tree);

  /* Duplicate labels are irrelevant to node-indexed C observations. */
  OK(crpx_tree_from_newick("(same:1,same:2);", &tree, NULL));
  size_t states[] = {SIZE_MAX, 0, 1};
  bt_parsimony_result *parsimony;
  OK(bt_parsimony_asr(tree, states, 2, false, &parsimony));
  close_value(parsimony->score, 1);
  OK(bt_discrete_asr(tree, states, 2, rate, false, &result));
  assert(!strcmp(crpx_tree_label(tree, 1), "same") && !strcmp(crpx_tree_label(tree, 2), "same"));
  crpx_tree_free(tree);
  close_value(parsimony->score, 1);
  close_value(result->probabilities[0] + result->probabilities[1], 1);
  bt_parsimony_free(parsimony); bt_discrete_free(result);
}

int main(void)
{
  exhaustive("((A:0.1,B:0.2):0.3,C:0.4);",2,false);
  exhaustive("((A:0.1,B:0.2):0.3,C:0.4);",3,true);
  exhaustive("(A:0.2,B:0.7,C:0.1);",4,false);
  exhaustive("(A:0,B:0.7,C:0.1);",3,true);
  exhaustive("((A:0.05):0.2,B:0.7,C:0.1);",3,true);
  special_cases(); underflow_case(); uncertain_observations();
  puts("discrete enumeration, uncertain observations, prepared rate fitting, parsimony and underflow: passed");
  return 0;
}
