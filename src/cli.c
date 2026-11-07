/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "obuina.h"
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { char *text; char ***row; size_t rows, columns; } table;
typedef struct { const char *name; size_t row; } named_row;
typedef struct {
  const char *command, *tree, *data, *trait, *predictor, *prefix, *traits;
  size_t states;
  double rate;
  bool fit_rate, ordered, midpoint, extend, allow_missing;
} options;

static void usage(FILE *stream)
{
  fputs("obuina " BT_VERSION " — ASR and comparative methods for labelled traits\n"
    "Usage: obuina COMMAND [OPTIONS]\n"
    "Commands: inspect, prepare, upgma, brownian, discrete, parsimony, pic, pgls\n"
    "  --tree FILE       Newick tree (required except upgma)\n"
    "  --data FILE       Tab-separated metadata, header includes taxon\n"
    "  --trait COLUMN    Trait column (default: value)\n"
    "  --traits C1,C2    Independent trait analyses; requires --prefix\n"
    "  --predictors C1,C2 Regression predictors or extra UPGMA features\n"
    "  --predictor COL   Backwards-compatible single predictor\n"
    "  --missing MODE    error (default) or estimate (allow missing observations)\n"
    "  --states K        Integer state codes 0..K-1, including unobserved states\n"
    "  --rate NUMBER     Equal-rates diffusion multiplier (default: 1)\n"
    "  --fit-rate        Fit the equal-rates CTMC multiplier within fixed bounds\n"
    "  --ordered         Ordered parsimony; requires --states\n"
    "  --midpoint        Explicitly midpoint-root before mapping tip metadata\n"
    "  --extend-tips     Explicitly extend terminal edges to equal tip depths\n"
    "  --prefix PATH     Write .nodes.tsv/.summary.tsv/.tree.nwk and method tables\n"
    "Without --prefix, results go to stdout and summary metrics to stderr.\n"
    "Discrete classes are inferred lexicographically unless --states is supplied.\n"
    "Unknown tokens: ?, NA or empty. Discrete unknowns are always permitted.\n"
    "With --missing estimate, Brownian predicts missing tips; PIC uses observed\n"
    "tips and PGLS uses complete rows. No filled values are treated as observations.\n"
    "Newick syntax does not establish biological rooting: supply a justified root.\n", stream);
}

static int number(const char *text, double *out)
{
  char *end;
  errno = 0;
  double value = strtod(text, &end);
  if (!*text || *end || errno || !isfinite(value)) return 0;
  *out = value;
  return 1;
}

static int integer(const char *text, size_t *out)
{
  char *end;
  if (!*text || *text == '-' || *text == '+') return 0;
  errno = 0;
  unsigned long long value = strtoull(text, &end, 10);
  if (*end || errno || value > SIZE_MAX) return 0;
  *out = (size_t)value;
  return 1;
}

static char *read_text(const char *path)
{
  FILE *stream = fopen(path, "rb");
  if (!stream) { perror(path); return NULL; }
  size_t used = 0, capacity = 4096;
  char *text = malloc(capacity);
  if (!text) { fclose(stream); return NULL; }
  for (;;) {
    if (used == capacity - 1) {
      if (capacity > SIZE_MAX / 2) { free(text); fclose(stream); return NULL; }
      capacity *= 2;
      char *next = realloc(text, capacity);
      if (!next) { free(text); fclose(stream); return NULL; }
      text = next;
    }
    size_t count = fread(text + used, 1, capacity - 1 - used, stream);
    if (memchr(text + used, 0, count)) {
      fprintf(stderr, "%s: embedded NUL byte\n", path);
      free(text); fclose(stream); return NULL;
    }
    used += count;
    if (!count) break;
  }
  if (ferror(stream)) { perror(path); free(text); fclose(stream); return NULL; }
  fclose(stream);
  text[used] = 0;
  return text;
}

static void table_free(table *data)
{
  if (data->row) for (size_t i = 0; i < data->rows; ++i) free(data->row[i]);
  free(data->row); free(data->text);
  memset(data, 0, sizeof(*data));
}

static int table_read(const char *path, table *data)
{
  data->text = read_text(path);
  if (!data->text || !*data->text) return 0;
  size_t lines = 1;
  for (char *p = data->text; *p; ++p) if (*p == '\n') ++lines;
  data->row = calloc(lines, sizeof(*data->row));
  if (!data->row) return 0;
  char *line = data->text;
  while (*line) {
    char *next = strchr(line, '\n');
    if (next) *next++ = 0;
    size_t length = strlen(line);
    if (length && line[length - 1] == '\r') line[--length] = 0;
    if (!length) { fprintf(stderr, "Blank metadata row\n"); return 0; }
    size_t fields = 1;
    for (char *p = line; *p; ++p) if (*p == '\t') ++fields;
    if (!data->rows) data->columns = fields;
    if (fields != data->columns) { fprintf(stderr, "Metadata row %zu has %zu fields; expected %zu\n",
      data->rows + 1, fields, data->columns); return 0; }
    char **row = calloc(fields, sizeof(*row));
    if (!row) return 0;
    data->row[data->rows++] = row;
    row[0] = line;
    size_t j = 1;
    for (char *p = line; *p; ++p) if (*p == '\t') { *p = 0; row[j++] = p + 1; }
    if (!next) break;
    line = next;
  }
  if (data->rows < 2) { fprintf(stderr, "Metadata needs a header and at least one taxon\n"); return 0; }
  for (size_t j = 0; j < data->columns; ++j)
    for (size_t k = 0; k < j; ++k) if (!strcmp(data->row[0][j], data->row[0][k])) {
      fprintf(stderr, "Duplicate metadata column: %s\n", data->row[0][j]); return 0;
    }
  return 1;
}

static size_t column(const table *data, const char *name)
{
  for (size_t j = 0; j < data->columns; ++j) if (!strcmp(name, data->row[0][j])) return j;
  fprintf(stderr, "Missing metadata column: %s\n", name);
  return SIZE_MAX;
}

static int unknown(const char *value)
{ return !*value || !strcmp(value, "?") || !strcmp(value, "NA"); }

static int compare_strings(const void *a, const void *b)
{ return strcmp(*(char *const *)a, *(char *const *)b); }
static int compare_rows(const void *a, const void *b)
{ return strcmp(((const named_row *)a)->name, ((const named_row *)b)->name); }

typedef struct { char *storage; char **items; size_t count; } column_list;
static void list_free(column_list *list) { free(list->items); free(list->storage); }
static int list_parse(const char *text, column_list *list)
{
  memset(list, 0, sizeof(*list));
  if (!text || !*text) return 0;
  size_t count = 1, length = strlen(text);
  for (const char *p = text; *p; ++p) if (*p == ',') ++count;
  list->storage = malloc(length + 1);
  list->items = calloc(count, sizeof(*list->items));
  if (!list->storage || !list->items) { list_free(list); memset(list, 0, sizeof(*list)); return 0; }
  memcpy(list->storage, text, length + 1);
  char *next = list->storage;
  for (size_t i = 0; i < count; ++i) {
    char *value = next, *comma = strchr(value, ',');
    if (comma) { *comma = 0; next = comma + 1; }
    while (*value == ' ' || *value == '\t') ++value;
    char *end = value + strlen(value);
    while (end > value && (end[-1] == ' ' || end[-1] == '\t')) *--end = 0;
    if (!*value) goto invalid;
    for (size_t j = 0; j < i; ++j) if (!strcmp(list->items[j], value)) goto invalid;
    list->items[list->count++] = value;
  }
  return 1;
invalid:
  fprintf(stderr, "Column lists must have nonempty, distinct names\n");
  list_free(list); memset(list, 0, sizeof(*list)); return 0;
}

static FILE *output_open(const char *prefix, const char *suffix, FILE *fallback)
{
  if (!prefix) return fallback;
  size_t a = strlen(prefix), b = strlen(suffix);
  if (a > SIZE_MAX - b - 1) return NULL;
  char *path = malloc(a + b + 1);
  if (!path) return NULL;
  memcpy(path, prefix, a); memcpy(path + a, suffix, b + 1);
  FILE *stream = fopen(path, "w");
  if (!stream) perror(path);
  free(path);
  return stream;
}

static int output_close(FILE *stream, FILE *fallback)
{
  if (stream == fallback) return fflush(stream) == 0;
  return fclose(stream) == 0;
}

static int write_tree(const crpx_tree *tree, const char *prefix, bool only)
{
  char *newick = NULL;
  if (crpx_tree_to_newick(tree, &newick) != CRPX_OK) return 0;
  FILE *stream = output_open(prefix, ".tree.nwk", stdout);
  int ok = stream && fprintf(stream, "%s\n", newick) > 0;
  crpx_tree_string_free(newick);
  if (stream) ok = output_close(stream, stdout) && ok;
  (void)only;
  return ok;
}

static int write_nodes(const crpx_tree *tree, const options *opt, const bt_brownian_result *brownian,
                        const bt_discrete_result *discrete, const bt_parsimony_result *parsimony,
                        char *const *labels, size_t states, const unsigned char *observed,
                        const unsigned char *used)
{
  FILE *stream = output_open(opt->prefix, ".nodes.tsv", stdout);
  if (!stream) return 0;
  fputs("node_id\tparent_id\tlabel\tis_tip\tbranch_length\tobserved\tused_for_fit"
    "\testimate\tvariance\tstate\toptimal_states", stream);
  for (size_t s = 0; discrete && s < states; ++s) fprintf(stream, "\tp_%zu", s);
  fputc('\n', stream);
  for (size_t i = 0; i < crpx_tree_node_count(tree); ++i) {
    crpx_node_id node = (crpx_node_id)i, parent = crpx_tree_parent(tree, node);
    const char *label = crpx_tree_label(tree, node);
    fprintf(stream, "%zu\t", i);
    if (parent != CRPX_NODE_NONE) fprintf(stream, "%u", parent);
    fprintf(stream, "\t%s\t%d\t", label ? label : "", crpx_tree_first_child(tree, node) == CRPX_NODE_NONE);
    double edge;
    if (parent != CRPX_NODE_NONE && crpx_tree_branch_length(tree, node, &edge) == CRPX_OK)
      fprintf(stream, "%.17g", edge);
    fputc('\t', stream);
    bool tip = crpx_tree_first_child(tree, node) == CRPX_NODE_NONE;
    if (tip && observed) fprintf(stream, "%d", observed[i] != 0);
    fputc('\t', stream);
    if (tip && used) fprintf(stream, "%d", used[i] != 0);
    fputc('\t', stream);
    if (brownian) fprintf(stream, "%.17g\t%.17g\t\t", brownian->mean[i], brownian->variance[i]);
    else if (discrete) {
      size_t best = 0;
      for (size_t s = 1; s < states; ++s)
        if (discrete->probabilities[i * states + s] > discrete->probabilities[i * states + best]) best = s;
      fprintf(stream, "\t\t%s\t", labels[best]);
    } else if (parsimony) {
      fprintf(stream, "\t\t%s\t", labels[parsimony->assignment[i]]);
      bool first = true;
      for (size_t s = 0; s < states; ++s) if (parsimony->optimal[i * states + s]) {
        fprintf(stream, "%s%zu", first ? "" : ",", s); first = false;
      }
    } else fputs("\t\t\t", stream);
    for (size_t s = 0; discrete && s < states; ++s) fprintf(stream, "\t%.17g", discrete->probabilities[i * states + s]);
    fputc('\n', stream);
  }
  int ok = !ferror(stream);
  return output_close(stream, stdout) && ok;
}

static int run_cli(int argc, char **argv);
static int run_batch(int argc, char **argv, const options *opt)
{
  if (!opt->prefix || !strcmp(opt->command, "inspect") || !strcmp(opt->command, "prepare") ||
      !strcmp(opt->command, "upgma")) {
    fprintf(stderr, "--traits requires an analysis command and --prefix\n"); return 2;
  }
  column_list list;
  if (!list_parse(opt->traits, &list)) return 2;
  char **arguments = calloc((size_t)argc + 1, sizeof(*arguments));
  size_t capacity = strlen(opt->prefix) + 40;
  char *prefix = malloc(capacity);
  FILE *manifest = NULL;
  int result = 0;
  if (!arguments || !prefix) { result = 1; goto done; }
  manifest = output_open(opt->prefix, ".traits.tsv", stderr);
  if (!manifest) { result = 1; goto done; }
  fputs("trait\tprefix\tstatus\n", manifest);
  for (size_t j = 0; j < list.count; ++j) {
    snprintf(prefix, capacity, "%s.trait%zu", opt->prefix, j + 1);
    for (int i = 0; i < argc; ++i) arguments[i] = argv[i];
    for (int i = 2; i + 1 < argc; ++i) {
      if (!strcmp(arguments[i], "--traits")) { arguments[i] = "--trait"; arguments[++i] = list.items[j]; }
      else if (!strcmp(arguments[i], "--prefix")) arguments[++i] = prefix;
    }
    int status = run_cli(argc, arguments);
    fprintf(manifest, "%s\t%s\t%s\n", list.items[j], prefix, status ? "error" : "ok");
    if (status) result = 1;
  }
  if (ferror(manifest)) result = 1;
done:
  if (manifest && !output_close(manifest, stderr)) result = 1;
  free(arguments); free(prefix); list_free(&list);
  return result;
}

static int run_cli(int argc, char **argv)
{
  options opt = {0};
  bool explicit_trait = false;
  opt.trait = "value"; opt.rate = 1;
  if (argc < 2 || !strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")) { usage(stdout); return argc < 2 ? 2 : 0; }
  if (!strcmp(argv[1], "--version")) { puts(bt_version()); return 0; }
  opt.command = argv[1];
  for (int i = 2; i < argc; ++i) {
    const char *arg = argv[i];
    if (!strcmp(arg, "--fit-rate")) opt.fit_rate = true;
    else if (!strcmp(arg, "--ordered")) opt.ordered = true;
    else if (!strcmp(arg, "--midpoint")) opt.midpoint = true;
    else if (!strcmp(arg, "--extend-tips")) opt.extend = true;
    else if (!strcmp(arg, "--help")) { usage(stdout); return 0; }
    else {
      if (++i >= argc) { fprintf(stderr, "Missing value for %s\n", arg); return 2; }
      const char *value = argv[i];
      if (!strcmp(arg, "--tree")) opt.tree = value;
      else if (!strcmp(arg, "--data")) opt.data = value;
      else if (!strcmp(arg, "--trait")) { opt.trait = value; explicit_trait = true; }
      else if (!strcmp(arg, "--traits")) opt.traits = value;
      else if (!strcmp(arg, "--predictor") || !strcmp(arg, "--predictors")) opt.predictor = value;
      else if (!strcmp(arg, "--prefix")) opt.prefix = value;
      else if (!strcmp(arg, "--missing") && !strcmp(value, "estimate")) opt.allow_missing = true;
      else if (!strcmp(arg, "--missing") && !strcmp(value, "error")) opt.allow_missing = false;
      else if (!strcmp(arg, "--states") && integer(value, &opt.states) && opt.states >= 1) {}
      else if (!strcmp(arg, "--rate") && number(value, &opt.rate) && opt.rate > 0) {}
      else { fprintf(stderr, "Invalid option or value: %s %s\n", arg, value); return 2; }
    }
  }
  if (opt.traits && explicit_trait) {
    fprintf(stderr, "Choose either --trait or --traits\n"); return 2;
  }
  if (opt.traits) return run_batch(argc, argv, &opt);
  bool infer = !strcmp(opt.command, "upgma"), prep = !strcmp(opt.command, "prepare");
  bool inspect = !strcmp(opt.command, "inspect"), bm = !strcmp(opt.command, "brownian");
  bool disc = !strcmp(opt.command, "discrete"), pars = !strcmp(opt.command, "parsimony");
  bool pic = !strcmp(opt.command, "pic"), gls = !strcmp(opt.command, "pgls");
  if (!(infer || prep || inspect || bm || disc || pars || pic || gls) || (!infer && !opt.tree)
      || (!(prep || inspect) && !opt.data) || (gls && !opt.predictor)
      || (opt.ordered && (!pars || !opt.states)) || (opt.fit_rate && !disc) || (disc && opt.states == 1)
      || (opt.states && !(disc || pars)) || (opt.predictor && !(gls || infer))) {
    usage(stderr); return 2;
  }
  crpx_tree *tree = NULL;
  table data = {0};
  double *values = NULL, *predictors = NULL, *metadata = NULL;
  size_t *observations = NULL;
  named_row *row_index = NULL;
  unsigned char *matched = NULL;
  unsigned char *observed_flags = NULL, *used_flags = NULL;
  column_list predictor_names = {0};
  size_t *predictor_columns = NULL;
  size_t observed_tips = 0, used_tips = 0;
  char **labels = NULL;
  size_t states = 0, allocated_labels = 0;
  bt_brownian_result *br = NULL;
  bt_discrete_result *dr = NULL;
  bt_parsimony_result *pr = NULL;
  bt_pic_result *ir = NULL;
  bt_pgls_result *gr = NULL;
  crpx_status status = CRPX_OK;
  int result = 1;
  double added = 0;
  size_t taxon_col = SIZE_MAX, trait_col = SIZE_MAX;
  if (opt.data) {
    if (!table_read(opt.data, &data)) goto cleanup;
    taxon_col = column(&data, "taxon"); trait_col = column(&data, opt.trait);
    if (taxon_col == SIZE_MAX || trait_col == SIZE_MAX) goto cleanup;
    if (opt.predictor) {
      if (!list_parse(opt.predictor, &predictor_names)) goto cleanup;
      predictor_columns = calloc(predictor_names.count, sizeof(*predictor_columns));
      if (!predictor_columns) { status = CRPX_ERROR_ALLOCATION; goto cleanup; }
      for (size_t j = 0; j < predictor_names.count; ++j) {
        predictor_columns[j] = column(&data, predictor_names.items[j]);
        if (predictor_columns[j] == SIZE_MAX) goto cleanup;
      }
    }
    row_index = calloc(data.rows - 1, sizeof(*row_index));
    matched = calloc(data.rows, sizeof(*matched));
    if (!row_index || !matched) { status = CRPX_ERROR_ALLOCATION; goto cleanup; }
    for (size_t i = 1; i < data.rows; ++i) {
      if (!*data.row[i][taxon_col]) { fprintf(stderr, "Empty metadata taxon\n"); goto cleanup; }
      row_index[i - 1] = (named_row){ data.row[i][taxon_col], i };
    }
    qsort(row_index, data.rows - 1, sizeof(*row_index), compare_rows);
    for (size_t i = 1; i < data.rows - 1; ++i)
      if (!strcmp(row_index[i - 1].name, row_index[i].name)) {
        fprintf(stderr, "Duplicate metadata taxon: %s\n", row_index[i].name); goto cleanup;
      }
  }
  if (infer) {
    size_t features = 1 + predictor_names.count, taxa = data.rows - 1;
    if (features > SIZE_MAX / taxa / sizeof(double)) { status = CRPX_ERROR_OVERFLOW; goto cleanup; }
    const char **names = calloc(taxa, sizeof(*names));
    metadata = calloc(taxa, features * sizeof(*metadata));
    if (!names || !metadata) { free(names); status = CRPX_ERROR_ALLOCATION; goto cleanup; }
    for (size_t i = 0; i < taxa; ++i) {
      names[i] = data.row[i + 1][taxon_col];
      bool valid = number(data.row[i + 1][trait_col], &metadata[i * features]);
      for (size_t j = 0; j < predictor_names.count; ++j)
        valid = number(data.row[i + 1][predictor_columns[j]], &metadata[i * features + j + 1]) && valid;
      if (!valid) {
        fprintf(stderr, "UPGMA metadata must be finite numbers\n"); free(names); goto cleanup;
      }
    }
    status = bt_tree_upgma(names, taxa, features, metadata, &tree);
    free(names);
  } else {
    char *text = read_text(opt.tree);
    if (!text) goto cleanup;
    size_t offset = 0;
    status = crpx_tree_from_newick(text, &tree, &offset);
    free(text);
    if (status != CRPX_OK) { fprintf(stderr, "Newick error at byte %zu\n", offset); goto cleanup; }
  }
  if (status != CRPX_OK) goto cleanup;
  if (opt.midpoint) {
    crpx_tree *next = NULL;
    status = bt_tree_midpoint_root(tree, &next);
    if (status != CRPX_OK) goto cleanup;
    crpx_tree_free(tree); tree = next;
  }
  if (opt.extend) {
    crpx_tree *next = NULL;
    status = bt_tree_extend_tips(tree, &next, &added);
    if (status != CRPX_OK) goto cleanup;
    crpx_tree_free(tree); tree = next;
  }
  bt_tree_info info;
  status = bt_tree_inspect(tree, 1e-9, &info);
  if (pars && (status == CRPX_ERROR_NOT_FOUND || status == CRPX_ERROR_INVALID_ARGUMENT)) {
    status = crpx_tree_validate(tree);
    info = (bt_tree_info){ crpx_tree_node_count(tree), crpx_tree_leaf_count(tree), NAN, NAN, false };
  }
  if (status != CRPX_OK) goto cleanup;
  for (size_t i = 0; i < info.nodes; ++i) {
    const char *label = crpx_tree_label(tree, (crpx_node_id)i);
    if (label && strpbrk(label, "\t\r\n")) {
      fprintf(stderr, "TSV output requires labels without tabs or newlines\n"); goto cleanup;
    }
  }
  if (prep || infer) {
    if (!write_tree(tree, opt.prefix, true)) goto cleanup;
  } else if (!(inspect)) {
    values = calloc(info.nodes, sizeof(*values));
    if (predictor_names.count > SIZE_MAX / info.nodes / sizeof(*predictors)) {
      status = CRPX_ERROR_OVERFLOW; goto cleanup;
    }
    predictors = calloc(info.nodes, (predictor_names.count ? predictor_names.count : 1) * sizeof(*predictors));
    observations = malloc(info.nodes * sizeof(*observations));
    observed_flags = calloc(info.nodes, 1); used_flags = calloc(info.nodes, 1);
    if (!values || !predictors || !observations || !observed_flags || !used_flags) {
      status = CRPX_ERROR_ALLOCATION; goto cleanup;
    }
    if (!opt.allow_missing && data.rows - 1 != info.tips) {
      fprintf(stderr, "Metadata and tree need exactly the same taxa; use --missing estimate for omitted tips\n");
      goto cleanup;
    }
    for (size_t i = 0; i < info.nodes; ++i) {
      observations[i] = SIZE_MAX; values[i] = NAN;
      for (size_t j = 0; j < predictor_names.count; ++j) predictors[i * predictor_names.count + j] = NAN;
    }
    if (disc || pars) {
      if (opt.states) {
        states = opt.states;
        labels = calloc(states, sizeof(*labels));
        if (!labels) { status = CRPX_ERROR_ALLOCATION; goto cleanup; }
        for (size_t s = 0; s < states; ++s) {
          labels[s] = malloc(32);
          if (!labels[s]) { status = CRPX_ERROR_ALLOCATION; goto cleanup; }
          ++allocated_labels;
          snprintf(labels[s], 32, "%zu", s);
        }
      } else {
        labels = calloc(data.rows - 1, sizeof(*labels));
        if (!labels) { status = CRPX_ERROR_ALLOCATION; goto cleanup; }
        size_t observed = 0;
        for (size_t r = 1; r < data.rows; ++r) {
          char *value = data.row[r][trait_col];
          if (unknown(value)) continue;
          labels[observed++] = value;
        }
        qsort(labels, observed, sizeof(*labels), compare_strings);
        for (size_t s = 0; s < observed; ++s)
          if (!states || strcmp(labels[states - 1], labels[s])) labels[states++] = labels[s];
        if (!states || (disc && states < 2)) {
          fprintf(stderr, "Need observed classes (>=2 for CTMC), or explicit integer --states K\n"); goto cleanup;
        }
      }
    }
    for (size_t i = 0; i < info.nodes; ++i) {
      if (crpx_tree_first_child(tree, (crpx_node_id)i) != CRPX_NODE_NONE) continue;
      const char *name = crpx_tree_label(tree, (crpx_node_id)i);
      if (!name || !*name) { fprintf(stderr, "Named tree tips required\n"); goto cleanup; }
      named_row key = { name, 0 };
      const named_row *entry = bsearch(&key, row_index, data.rows - 1, sizeof(*row_index), compare_rows);
      if (!entry) {
        if (opt.allow_missing) continue;
        fprintf(stderr, "Missing metadata taxon: %s\n", name); goto cleanup;
      }
      size_t r = entry->row;
      if (matched[r]) { fprintf(stderr, "Duplicate tree tip: %s\n", name); goto cleanup; }
      matched[r] = 1;
      const char *value = data.row[r][trait_col];
      if (disc || pars) {
        if (unknown(value)) continue;
        if (opt.states) {
          if (!integer(value, &observations[i]) || observations[i] >= states) {
            fprintf(stderr, "State for %s must be in 0..%zu\n", name, states - 1); goto cleanup;
          }
        } else {
          char *key_label = (char *)value;
          char **entry_label = bsearch(&key_label, labels, states, sizeof(*labels), compare_strings);
          if (!entry_label) { status = CRPX_ERROR_INVALID_ARGUMENT; goto cleanup; }
          observations[i] = (size_t)(entry_label - labels);
        }
      } else if (unknown(value) && opt.allow_missing) {
        values[i] = NAN;
      } else if (!number(value, &values[i])) {
        fprintf(stderr, "Finite continuous observation required for %s; use --missing estimate for unknowns\n", name);
        goto cleanup;
      }
      observed_flags[i] = (disc || pars) ? observations[i] != SIZE_MAX : !isnan(values[i]);
      used_flags[i] = observed_flags[i];
      if (gls) for (size_t j = 0; j < predictor_names.count; ++j) {
        const char *input = data.row[r][predictor_columns[j]];
        if (unknown(input) && opt.allow_missing) used_flags[i] = 0;
        else if (!number(input, &predictors[i * predictor_names.count + j])) {
          fprintf(stderr, "Finite predictor required for %s\n", name); goto cleanup;
        }
      }
      observed_tips += observed_flags[i]; used_tips += used_flags[i];
    }
    for (size_t r = 1; r < data.rows; ++r) if (!matched[r]) {
      fprintf(stderr, "Metadata taxon absent from tree: %s\n", data.row[r][taxon_col]); goto cleanup;
    }
    if (bm) status = bt_brownian_fit(tree, values, &br);
    else if (disc) status = bt_discrete_asr(tree, observations, states, opt.rate, opt.fit_rate, &dr);
    else if (pars) status = bt_parsimony_asr(tree, observations, states, opt.ordered, &pr);
    else if (pic) status = bt_pic(tree, values, &ir);
    else if (gls) status = bt_pgls(tree, values, predictors, predictor_names.count, &gr);
    if (status != CRPX_OK) goto cleanup;
  }
  if (!(prep || infer) && !write_nodes(tree, &opt, br, dr, pr, labels, states, observed_flags, used_flags))
    goto cleanup;
  if (opt.prefix && !(prep || infer) && !write_tree(tree, opt.prefix, false)) goto cleanup;
  FILE *summary = output_open(opt.prefix, ".summary.tsv", stderr);
  if (!summary) goto cleanup;
  fprintf(summary, "metric\tvalue\nversion\t%s\nmethod\t%s\nnodes\t%zu\ntips\t%zu\nultrametric\t%s\n"
    "minimum_tip_depth\t%.17g\nmaximum_tip_depth\t%.17g\nmidpoint_requested\t%d\nterminal_extension\t%.17g\n",
    bt_version(), opt.command, info.nodes, info.tips,
    isfinite(info.maximum_tip_depth) ? (info.ultrametric ? "1" : "0") : "NA",
    info.minimum_tip_depth, info.maximum_tip_depth, opt.midpoint, added);
  if (br) fprintf(summary, "root_mean\t%.17g\nsigma2_ml\t%.17g\nlog_likelihood\t%.17g\n", br->root_mean, br->sigma2, br->log_likelihood);
  if (observed_flags) fprintf(summary, "trait\t%s\nmissing_policy\t%s\nobserved_tips\t%zu\nmissing_tips\t%zu\n"
    "used_tips\t%zu\n", opt.trait, opt.allow_missing ? "estimate" : "error", observed_tips,
    info.tips - observed_tips, used_tips);
  if (dr) fprintf(summary, "rate\t%.17g\nrate_fitted\t%d\nlog_likelihood\t%.17g\n", dr->rate, opt.fit_rate, dr->log_likelihood);
  if (pr) fprintf(summary, "parsimony_score\t%.17g\nordered\t%d\n", pr->score, opt.ordered);
  if (gr) fprintf(summary, "sigma2_ml\t%.17g\nlog_likelihood\t%.17g\n", gr->sigma2, gr->log_likelihood);
  if (ir) fprintf(summary, "contrasts\t%zu\n", ir->count);
  for (size_t s = 0; s < states; ++s) fprintf(summary, "state_%zu\t%s\n", s, labels[s]);
  int ok = !ferror(summary);
  if (!output_close(summary, stderr) || !ok) goto cleanup;
  if (ir) {
    FILE *stream = output_open(opt.prefix, ".contrasts.tsv", stdout);
    if (!stream) goto cleanup;
    fputs("node_id\tstandardised_contrast\tcontrast_variance\n", stream);
    for (size_t i = 0; i < ir->count; ++i)
      fprintf(stream, "%u\t%.17g\t%.17g\n", ir->node[i], ir->contrast[i], ir->variance[i]);
    ok = !ferror(stream);
    if (!output_close(stream, stdout) || !ok) goto cleanup;
  }
  if (gr) {
    FILE *stream = output_open(opt.prefix, ".coefficients.tsv", stdout);
    if (!stream) goto cleanup;
    fputs("parameter\tcoefficient\tstandard_error\n", stream);
    for (size_t j = 0; j < gr->parameters; ++j)
      fprintf(stream, "%s\t%.17g\t%.17g\n", j ? predictor_names.items[j - 1] : "intercept",
        gr->coefficients[j], gr->standard_errors[j]);
    ok = !ferror(stream);
    if (!output_close(stream, stdout) || !ok) goto cleanup;
  }
  result = 0;
cleanup:
  if (status != CRPX_OK) fprintf(stderr, "obuina: %s (see docs/methods.md for method restrictions)\n", crpx_status_string(status));
  bt_brownian_free(br); bt_discrete_free(dr); bt_parsimony_free(pr); bt_pic_free(ir); bt_pgls_free(gr);
  crpx_tree_free(tree); table_free(&data);
  free(values); free(predictors); free(observations); free(metadata);
  free(row_index); free(matched);
  free(observed_flags); free(used_flags); free(predictor_columns); list_free(&predictor_names);
  for (size_t s = 0; s < allocated_labels; ++s) free(labels[s]);
  free(labels);
  return result;
}

int main(int argc, char **argv) { return run_cli(argc, argv); }
