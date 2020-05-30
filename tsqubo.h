#pragma once

#include <stdlib.h>
#include <string.h>

/** A QUBO solution structure. */
struct tsqubo_solution {
  /** The objective function value. */
  double fx;
  /** The solution vector. */
  double *x;
  /** The reevaluation vector. */
  double *dx;
  /** Permutation vector which contains the indices 1,2,...,n. */
  size_t *N;
  /** Associates each variable to a position of @ref N. */
  size_t *I;
};

/**
 * Initializes the TS structure.
 *
 * @param sol an uninitialized solution structure.
 * @param n the number of variables.
 */
void tsqubo_solution_init(struct tsqubo_solution *sol, size_t n) {
  sol->fx = 0;
  sol->x = (double *)calloc(n, sizeof(double));
  sol->dx = (double *)calloc(n, sizeof(double));
  sol->N = (size_t *)calloc(n, sizeof(size_t));
  sol->I = (size_t *)calloc(n, sizeof(size_t));
  for (size_t i = 0; i < n; i++) sol->N[i] = sol->I[i] = i;
}

/**
 * Frees the TS structure.
 *
 * @param sol an initialized solution structure.
 */
void tsqubo_solution_free(struct tsqubo_solution *sol) {
  free(sol->x);
  free(sol->dx);
  free(sol->N);
  free(sol->I);
}

/**
 * Assigns a position, in the permuation vector of a solution, to the given variable.
 *
 * @param sol an initialized solution structure.
 * @param i the index of the variable.
 * @param pos the assigned position in the permutation vector.
 */
void tsqubo_solution_reorder(struct tsqubo_solution *sol, size_t i, size_t pos) {
  sol->I[sol->N[pos]] = sol->I[i];
  sol->N[sol->I[i]] = sol->N[pos];
  sol->I[i] = pos;
  sol->N[pos] = i;
}

/**
 * Searches for the best flip within the best first `nvars` variables
 * in the permutation vector of a solution.
 *
 * @param sol an initialized solution structure.
 * @param nvars the number of variables to search, which must be greater than 0.
 * @return the index of the chosen variable to flip.
 */
size_t tsqubo_solution_find_flip(const struct tsqubo_solution *sol, size_t nvars) {
  size_t i = sol->N[0];
  double d = sol->dx[i];
  for (size_t k = 1; k < nvars; k++) {
    size_t j = sol->N[k];
    double e = sol->dx[j];
    if (e < d) {
      d = e;
      i = j;
    }
  }
  return i;
}

/** A QUBO tabu list structure. */
struct tsqubo_tabulist {
  /** The tabu list. */
  size_t *L;
  /** Specifies that the first @ref nnontabu indices are not tabu. */
  size_t nnontabu;
  /**
   * Specifies that the first @ref ntargets indices either
   * are not tabu or meet the aspiration criterion.
   */
  size_t ntargets;
};

/**
 * Initializes the tabu list structure.
 *
 * @param tabu an uninitialized tabu list structure.
 * @param n the number of variables in the problem instance.
 */
void tsqubo_tabulist_init(struct tsqubo_tabulist *tabu, size_t n) {
  tabu->L = (size_t *)calloc(n, sizeof(size_t));
  tabu->nnontabu = n;
  tabu->ntargets = n;
}

/**
 * Frees memory associated with a tabu list structure.
 *
 * @param tabu an initialized tabu list structure.
 */
void tsqubo_tabulist_free(struct tsqubo_tabulist *tabu) { free(tabu->L); }

/** A component of the Q matrix of a QUBO problem instance. */
struct tsqubo_instance_component {
  /** The row of the component. */
  size_t i;
  /** The column of the component. */
  size_t j;
  /** The value of the component. */
  double q;
};

/** A QUBO problem instance structure. */
struct tsqubo_instance {
  /** The current number of variables in the instance. */
  size_t n;
  /** The array of components in the matrix. */
  struct tsqubo_instance_component *components;
  /** The current size of @ref components. */
  size_t ncomponents;
  /** The capacity of @ref components. */
  size_t maxcomponents;
};

/**
 * Initializes the problem instance structure.
 *
 * @param inst an uninitialized problem instance structure.
 * @param capacity the initial capacity for components in the matrix.
 */
void tsqubo_instance_init(struct tsqubo_instance *inst, size_t capacity) {
  inst->n = 0;
  inst->components = (struct tsqubo_instance_component *)malloc(
      (inst->maxcomponents = (capacity < 2 ? 2 : capacity)) *
      sizeof(struct tsqubo_instance_component));
  inst->ncomponents = 0;
}

/**
 * Allocates and initializes a problem instance structure.
 *
 * @param capacity the initial capacity for components in the matrix.
 * @return a newly-allocated QUBO problem instance structure.
 */
struct tsqubo_instance *tsqubo_instance_new(size_t capacity) {
  struct tsqubo_instance *inst =
      (struct tsqubo_instance *)calloc(1, sizeof(struct tsqubo_instance));
  tsqubo_instance_init(inst, capacity);
  return inst;
}

/**
 * Frees the problem instance structure.
 *
 * @param inst an initialized problem instance structure.
 */
void tsqubo_instance_free(struct tsqubo_instance *inst) { free(inst->components); }

/**
 * Sets a coefficient of the Q matrix.
 *
 * The coefficient should NOT have been previously set,
 * otherwise performance can be affected.
 *
 * @param inst an initialized problem instance structure.
 * @param i the row of the matrix.
 * @param j the column of the matrix.
 * @param q the value of the coefficient.
 */
void tsqubo_instance_add_component(struct tsqubo_instance *inst, size_t i, size_t j, double q) {
  // adjust n if needed
  if (i >= inst->n)
    inst->n = i + 1;
  else if (j >= inst->n)
    inst->n = j + 1;
  // resize array of components if needed
  if (inst->ncomponents + 2 > inst->maxcomponents)
    inst->components = (struct tsqubo_instance_component *)realloc(
        inst->components, (inst->maxcomponents *= 2) * sizeof(struct tsqubo_instance_component));
  // assign values to component
  if (i == j) {
    // diagonal component
    inst->components[inst->ncomponents].i = i;
    inst->components[inst->ncomponents].j = j;
    inst->components[inst->ncomponents].q = q;
    inst->ncomponents++;
  } else {
    // non-diagonal component
    inst->components[inst->ncomponents].i = i;
    inst->components[inst->ncomponents].j = j;
    inst->components[inst->ncomponents].q = 2 * q;
    inst->ncomponents++;
    inst->components[inst->ncomponents].i = j;
    inst->components[inst->ncomponents].j = i;
    inst->components[inst->ncomponents].q = 2 * q;
    inst->ncomponents++;
  }
}

/**
 * Compares two components of a QUBO problem instance.
 *
 * A diagonal component should always come before non-diagonal components.
 */
static int compare_instance_components(const struct tsqubo_instance_component *a,
                                       const struct tsqubo_instance_component *b) {
  if (a->i < b->i) return -1;
  if (b->i < a->i) return 1;
  if (a->j == a->i) return -1;
  if (b->j == b->i) return 1;
  if (a->j < b->j) return -1;
  if (b->j < a->j) return 1;
  return 0;
}

/** A QUBO problem instance structure, optimized for quick access. */
struct tsqubo_compiled_instance {
  /** The number of variables, or rows in @ref Q. */
  size_t n;
  /** The coefficients of the matrix. */
  double *Q;
  /** The columns of rows. */
  size_t *cols;
  /** The position in @ref cols containing the columns of each row. */
  size_t *icols;
  /** The number of columns of each row in @ref cols after @ref icols. */
  size_t *ncols;
};

void tsqubo_compiled_instance_init(struct tsqubo_compiled_instance *cinst,
                                   struct tsqubo_instance *inst) {
  qsort(inst->components, inst->ncomponents, sizeof(struct tsqubo_instance_component),
        (int (*)(const void *, const void *))compare_instance_components);
  cinst->n = inst->n;
  cinst->Q = (double *)calloc(inst->ncomponents, sizeof(double));
  cinst->cols = (size_t *)calloc(inst->ncomponents, sizeof(size_t));
  cinst->icols = (size_t *)calloc(inst->n, sizeof(size_t));
  cinst->ncols = (size_t *)calloc(inst->n, sizeof(size_t));
  for (size_t k = 0; k < inst->ncomponents; k++) {
    size_t i = inst->components[k].i, j = inst->components[k].j;
    cinst->Q[k] = inst->components[k].q;
    cinst->cols[k] = j;
    if (k > 0 && i != inst->components[k - 1].i) cinst->icols[i] = k;
    cinst->ncols[i]++;
  }
}

struct tsqubo_compiled_instance *tsqubo_compiled_instance_new(struct tsqubo_instance *inst) {
  struct tsqubo_compiled_instance *cinst =
      (struct tsqubo_compiled_instance *)calloc(1, sizeof(struct tsqubo_compiled_instance));
  tsqubo_compiled_instance_init(cinst, inst);
  return cinst;
}

void tsqubo_compiled_instance_free(struct tsqubo_compiled_instance *cinst) {
  free(cinst->Q);
  free(cinst->cols);
  free(cinst->icols);
  free(cinst->ncols);
}

/**
 * Copies a solution.
 *
 * @param cinst the instance which the solutions belong to.
 * @param dst the destination solution.
 * @param src the source solution.
 */
void tsqubo_compiled_instance_copy_solution(const struct tsqubo_compiled_instance *cinst,
                                            struct tsqubo_solution *dst,
                                            const struct tsqubo_solution *src) {
  dst->fx = src->fx;
  memcpy(dst->x, src->x, cinst->n * sizeof(double));
  memcpy(dst->dx, src->dx, cinst->n * sizeof(double));
}

/**
 * Flips a variable of a solution.
 *
 * @param inst an initialized problem instance structure.
 * @param sol an initialized solution structure.
 * @param i the index of the variable to flip.
 */
void tsqubo_compiled_instance_flip_solution(const struct tsqubo_compiled_instance *inst,
                                            struct tsqubo_solution *sol, size_t i) {
  sol->fx += sol->dx[i];
  sol->x[i] = !sol->x[i];
  sol->dx[i] = -sol->dx[i];
  for (size_t k = 1; k < inst->ncols[i]; k++) {
    size_t j = inst->cols[inst->icols[i] + k];
    sol->dx[j] -= (1 - 2 * sol->x[i]) * (1 - 2 * sol->x[j]) * inst->Q[inst->icols[i] + k];
  }
}

/** A QUBO tabu search structure. */
struct tsqubo {
  /** The problem instance. */
  struct tsqubo_compiled_instance inst;
  /** The incumbent solution of the search. */
  struct tsqubo_solution inc;
  /** The current solution of the search. */
  struct tsqubo_solution cur;
  /** The tabu list. */
  struct tsqubo_tabulist tabu;
};

/**
 * Initializes the TS structure.
 *
 * @param ts an uninitialized TS structure.
 */
void tsqubo_init(struct tsqubo *ts, struct tsqubo_instance *inst) {
  tsqubo_compiled_instance_init(&ts->inst, inst);
  tsqubo_solution_init(&ts->inc, inst->n);
  tsqubo_solution_init(&ts->cur, inst->n);
  tsqubo_tabulist_init(&ts->tabu, inst->n);
}

/**
 * Allocates and initializes the TS structure.
 *
 * @return a newly-allocated QUBO problem instance structure.
 */
struct tsqubo *tsqubo_new(struct tsqubo_instance *inst) {
  struct tsqubo *ts = (struct tsqubo *)calloc(1, sizeof(struct tsqubo));
  tsqubo_init(ts, inst);
  return ts;
}

/**
 * Frees memory associated with a TS structure.
 *
 * @param ts an initialized TS structure.
 */
void tsqubo_free(struct tsqubo *ts) {
  tsqubo_compiled_instance_free(&ts->inst);
  tsqubo_solution_free(&ts->inc);
  tsqubo_solution_free(&ts->cur);
  tsqubo_tabulist_free(&ts->tabu);
}

/**
 * Copies the current solution to the incumbent solution.
 *
 * @param ts an initialized TS structure.
 */
void tsqubo_commit_incumbent(struct tsqubo *ts) {
  tsqubo_compiled_instance_copy_solution(&ts->inst, &ts->inc, &ts->cur);
}

/**
 * Zeroes the current and incumbent solution.
 *
 * @param ts an initialized TS structure.
 */
void tsqubo_reset_solutions(struct tsqubo *ts) {
  ts->cur.fx = 0;
  memset(ts->cur.x, 0, ts->inst.n * sizeof(double));
  for (size_t i = 0; i < ts->inst.n; i++) ts->cur.dx[i] = ts->inst.Q[ts->inst.icols[i]];
  tsqubo_commit_incumbent(ts);
}

/**
 * Resets the tabu list.
 *
 * @param ts an initialized TS structure.
 */
void tsqubo_reset_tabu(struct tsqubo *ts) {
  ts->tabu.nnontabu = ts->inst.n;
  ts->tabu.ntargets = ts->inst.n;
  memset(ts->tabu.L, 0, ts->inst.n * sizeof(size_t));
}

/**
 * Flips a variable of the current solution.
 *
 * @param ts an initialized TS structure.
 * @param i the index of the variable to flip.
 * @return whether the incumbent solution was improved.
 */
int tsqubo_flip_current(struct tsqubo *ts, size_t i) {
  tsqubo_compiled_instance_flip_solution(&ts->inst, &ts->cur, i);
  if (ts->cur.fx >= ts->inc.fx) return 0;
  tsqubo_commit_incumbent(ts);
  return 1;
}

/**
 * Performs a local search on the incumbent solution of the search.
 *
 * @param ts an initialized TS structure.
 */
void tsqubo_local_search(struct tsqubo *ts) {
  // partition the N vector between variables which offer improvement when flipped or not
  size_t nvars = 0;
  for (size_t k = 0; k < ts->inst.n; k++) {
    size_t i = ts->inc.N[k];
    if (ts->inc.dx[i] < 0) tsqubo_solution_reorder(&ts->inc, i, nvars++);
  }
  // local search
  while (nvars) {
    size_t i = tsqubo_solution_find_flip(&ts->inc, nvars);
    tsqubo_compiled_instance_flip_solution(&ts->inst, &ts->inc, i);
    tsqubo_solution_reorder(&ts->inc, i, --nvars);
    for (size_t k = 1; k < ts->inst.ncols[i]; k++) {
      size_t j = ts->inst.cols[ts->inst.icols[i] + k];
      if (ts->inc.dx[j] < 0 && ts->inc.I[j] >= nvars)
        tsqubo_solution_reorder(&ts->inc, j, nvars++);
      else if (ts->inc.dx[j] >= 0 && ts->inc.I[j] < nvars)
        tsqubo_solution_reorder(&ts->inc, j, --nvars);
    }
  }
  // update in case some variables start or stop meeting the aspiration criterion
  size_t k = ts->tabu.nnontabu, ntargets_prev = ts->tabu.ntargets;
  while (k < ts->tabu.ntargets) {
    size_t j = ts->cur.N[k];
    if (ts->cur.fx + ts->cur.dx[j] >= ts->inc.fx)
      tsqubo_solution_reorder(&ts->cur, j, --ts->tabu.ntargets);
    else
      k++;
  }
  for (k = ntargets_prev; k < ts->inst.n; k++) {
    size_t j = ts->cur.N[k];
    if (ts->cur.fx + ts->cur.dx[j] < ts->inc.fx)
      tsqubo_solution_reorder(&ts->cur, j, ts->tabu.ntargets++);
  }
}

/** The state of the tabu search. */
enum tsqubo_state {
  /** No flip moves could be made. */
  TSQUBO_STATE_END,
  /** A flip move was made but did not improve the incumbent. */
  TSQUBO_STATE_SEARCH,
  /** A flip move was made which improved the incumbent. */
  TSQUBO_STATE_IMPROVEMENT
};

/**
 * Performs a TS iteration.
 *
 * @param ts an initialized TS structure.
 * @param K the tabu tenure to assign to flipped variables.
 * @return the result.
 */
enum tsqubo_state tsqubo_iterate(struct tsqubo *ts, size_t K) {
  if (!ts->tabu.ntargets) return TSQUBO_STATE_END;
  size_t i = tsqubo_solution_find_flip(&ts->cur, ts->tabu.ntargets);
  int improvement = tsqubo_flip_current(ts, i);
  // decrement aspiration and tabu variables
  for (size_t k = ts->tabu.nnontabu; k < ts->inst.n; k++) {
    size_t j = ts->cur.N[k], l = --ts->tabu.L[j];
    if (k >= ts->tabu.ntargets && (!l || ts->cur.fx + ts->cur.dx[j] < ts->inc.fx))
      tsqubo_solution_reorder(&ts->cur, j, ts->tabu.ntargets++);
    if (!l) tsqubo_solution_reorder(&ts->cur, j, ts->tabu.nnontabu++);
  }
  // assign tenure to flipped variable
  ts->tabu.L[i] = K;
  if (ts->cur.I[i] < ts->tabu.nnontabu) tsqubo_solution_reorder(&ts->cur, i, --ts->tabu.nnontabu);
  if (ts->cur.fx + ts->cur.dx[i] >= ts->inc.fx)
    tsqubo_solution_reorder(&ts->cur, i, --ts->tabu.ntargets);
  // update variables whose delta values were updated
  for (size_t k = 1; k < ts->inst.ncols[i]; k++) {
    size_t j = ts->inst.cols[ts->inst.icols[i] + k];
    if (ts->cur.I[j] >= ts->tabu.ntargets) {
      if (ts->cur.fx + ts->cur.dx[j] < ts->inc.fx)
        tsqubo_solution_reorder(&ts->cur, j, ts->tabu.ntargets++);
    } else if (ts->cur.I[j] >= ts->tabu.nnontabu) {
      if (ts->cur.fx + ts->cur.dx[j] >= ts->inc.fx)
        tsqubo_solution_reorder(&ts->cur, j, --ts->tabu.ntargets);
    }
  }
  return improvement ? TSQUBO_STATE_IMPROVEMENT : TSQUBO_STATE_SEARCH;
}

/**
 * Performs TS iterations until no improvements are possible for a number of iterations.
 *
 * @param ts an initialized TS structure.
 * @param K the tabu tenure to assign to flipped variables.
 * @param cutoff the improvement cutoff.
 */
void tsqubo_iterate_cutoff(struct tsqubo *ts, size_t K, size_t cutoff) {
begin:
  for (size_t i = 0; i < cutoff; i++) {
    switch (tsqubo_iterate(ts, K)) {
      case TSQUBO_STATE_IMPROVEMENT:
        tsqubo_local_search(ts);
        goto begin;
      case TSQUBO_STATE_END:
        return;
      default:
        break;
    }
  }
}