# Curupixa findings and adoption status

Obuina 0.2.0 adopts curupixa commit **`a1adfc4258fe9c6bc771464934db9d2ef6d5d397`**, which follows the initial findings recorded during Obuina 0.1.0 development against `268040a6f359b7b3d98f346211076559d667c6b6`. The upstream [downstream-trait report](../../curupixa/docs/downstream-trait-api.md) describes the changes and their validation. These findings concern public API coverage, modelling policy and performance; no confirmed numerical correctness defect was identified in the original dependency inspection.

## 1. Rooting inside an edge: upstream primitive added

Curupixa now exposes `crpx_tree_reroot_edge`, with explicit policies for reusing an artificial binary root or preserving the old vertex and appending a root. It preserves existing IDs and labels, supports endpoint splits and documents unary-root restrictions. See [tree extensions](../../curupixa/docs/tree-extensions.md) and [tree.h](../../curupixa/lib/tree.h).

Obuina retains its midpoint diameter search and graph reconstruction because its established transformation suppresses unlabelled degree-two vertices and removes unary root stems. The upstream primitive deliberately preserves vertices and rejects unary-root chains; substituting it directly would change labels, vertex counts and the existing midpoint contract. Node IDs can therefore still change after `bt_tree_midpoint_root`. A future explicit edge/outgroup-rooting application API can use the new primitive without imposing midpoint's contraction policy.

## 2. Broad-polytomy appends: builder adopted

The ordinary upstream `crpx_tree_add_node` remains convenient for small bounded-degree construction. Its sibling scanning is intentional. The new optional `crpx_tree_builder` caches child tails outside the compact tree and gives amortised constant-time ordered appends after preparation. Labels/lengths can be set between appends without invalidating that cache.

Obuina midpoint rebuilding now uses the builder for every appended node, removing its former quadratic sibling-scan bottleneck on star trees. Builder cleanup precedes destruction of the bound tree. Midpoint search, contraction and reconstruction are linear in node count; the complete operation also includes tip-name validation with $O(T\log T)$ sorting. A 4,096-tip annotated-star regression checks distances, midpoint depths, labels, node counts and unchanged input. No machine-specific timing threshold is imposed.

## 3. Full reconstruction and continuous kernels: upstream APIs adopted

Curupixa's [reconstruction API](../../curupixa/docs/reconstruction.md) now computes full-data marginal state probabilities and coherent parsimony assignments/global admissible-state sets. Its [continuous-trait API](../../curupixa/docs/continuous-traits.md) supplies linear-memory Brownian ML, optional ancestral marginals and PIC with explicit tip-node mappings.

Obuina delegates those numerical kernels through its existing C result wrappers. It retains table/column handling, batch interfaces, application rate-search policy, PGLS and missing-continuous-observation handling. The upstream continuous routines still require complete observations for the tree passed to them; Obuina builds the appropriate observed subtree and maps the resulting estimates/predictions back to original IDs. This keeps the original-root likelihood convention while adding downstream missing-tip support.

## 4. Repeated likelihood setup: prepared binding adopted

Curupixa now provides `crpx_likelihood_binding_new` and `crpx_likelihood_evaluate_prepared`. A binding owns a tree snapshot and retains immutable observation data. After warm-up, fixed-dimension rate-only evaluation avoids repeated mapping, structural validation, traversal construction and capacity allocations. The ordinary one-off interface remains available.

Obuina prepares a binding once for the equal-rates scalar search and requests reconstruction at the final rate. Reusable public Obuina workspaces, cross-trait batching with shared preparation and complete end-to-end allocation benchmarks remain future opportunities. Prepared snapshots have explicit lifetimes; changing an original tree does not silently update a previously prepared binding.

## 5. Operation-specific validation: shared policy adopted

`crpx_tree_validate_policy` now exposes named/unique-tip, complete/nonnegative/positive-length, binary-degree and root-stem checks. Obuina tree inspection uses this shared policy for its named unique tips and complete nonnegative edges, then calculates its depth and ultrametric summaries. The root stem remains excluded.

Curupixa's broad tree representation still permits missing/signed lengths and arbitrary labels because other algorithms use those representations. Obuina' probabilistic models apply their own requirements; parsimony ignores lengths. These are operation-specific modelling contracts rather than storage defects.

## Remaining downstream work

The [limitations document](limitations.md) tracks richer evolutionary models, dense PGLS scaling, continuous measurement-error likelihoods and joint multivariate traits. The [missing-data and uncertainty guide](missing-data-and-uncertainty.md) specifies what is now implemented and which observation-model/API choices remain proposals. The upstream source was not modified during this Obuina adoption.
