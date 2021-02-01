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
}

/**
 * Frees the TS structure.
 *
 * @param sol an initialized solution structure.
 */
void tsqubo_solution_free(struct tsqubo_solution *sol) {
  free(sol->x);
  free(sol->dx);
}

/** A QUBO problem instance structure. */
struct tsqubo_instance {
  /** The current number of variables in the instance. */
  size_t n;
  /** The array of components in the matrix. */
  struct {
    /** The row of the component. */
    size_t i;
    /** The column of the component. */
    size_t j;
    /** The value of the component. */
    double q;
  } * components;
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
  inst->components = (__typeof__(inst->components))malloc(
      (inst->maxcomponents = (capacity < 2 ? 2 : capacity)) * sizeof(inst->components[0]));
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
    inst->components = (__typeof__(inst->components))realloc(
        inst->components, (inst->maxcomponents *= 2) * sizeof(inst->components[0]));
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
static int compare_instance_components(const void *a_, const void *b_) {
  const __typeof__(((struct tsqubo_instance *)NULL)->components)
      a = (__typeof__(((struct tsqubo_instance *)NULL)->components))a_,
      b = (__typeof__(((struct tsqubo_instance *)NULL)->components))b_;
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
  size_t *C;
  /** The position in @ref cols containing the columns of each row. */
  size_t *R;
};

static void tsqubo_compiled_instance_init(struct tsqubo_compiled_instance *cinst,
                                          struct tsqubo_instance *inst) {
  qsort(inst->components, inst->ncomponents,
        sizeof(((struct tsqubo_instance *)NULL)->components[0]), compare_instance_components);
  cinst->n = inst->n;
  cinst->Q = (double *)calloc(inst->ncomponents, sizeof(double));
  cinst->C = (size_t *)calloc(inst->ncomponents, sizeof(size_t));
  cinst->R = (size_t *)calloc(inst->n + 1, sizeof(size_t));
  for (size_t k = 0; k < inst->ncomponents; k++) {
    size_t i = inst->components[k].i, j = inst->components[k].j;
    cinst->Q[k] = inst->components[k].q;
    cinst->C[k] = j;
    if (k > 0 && i != inst->components[k - 1].i) cinst->R[i] = k;
  }
  cinst->R[inst->n] = inst->ncomponents;
  return;
}

static void tsqubo_compiled_instance_free(struct tsqubo_compiled_instance *cinst) {
  free(cinst->Q);
  free(cinst->C);
  free(cinst->R);
}

/** A QUBO tabu search indexed priority queue structure. */
struct tsqubo_queue {
  /** A permutation of `{1,2,...n}`. */
  size_t *N;
  /** An index vector that associates each `{1,2,...n}` to its position in @ref N. */
  size_t *I;
  /** The size of the queue. */
  size_t size;
};

static void queue_init(struct tsqubo_queue *q, size_t n) {
  q->N = (size_t *)calloc(n, sizeof(size_t));
  q->I = (size_t *)calloc(n, sizeof(size_t));
  for (size_t i = 0; i < n; i++) q->N[i] = q->I[i] = i;
  q->size = 0;
}

static void queue_free(struct tsqubo_queue *q) {
  free(q->N);
  free(q->I);
}

static void queue_reorder(struct tsqubo_queue *q, size_t i, size_t k) {
  q->I[q->N[k]] = q->I[i];
  q->N[q->I[i]] = q->N[k];
  q->I[i] = k;
  q->N[k] = i;
}

static int queue_contains(const struct tsqubo_queue *q, size_t i) { return q->I[i] < q->size; }

static size_t queue_top(const struct tsqubo_queue *q) { return q->N[0]; }

static int compare_size(const struct tsqubo_queue *q, const void *v_, size_t a, size_t b) {
  const size_t *v = (const size_t *)v_;
  return v[q->N[a]] < v[q->N[b]] || (v[q->N[a]] == v[q->N[b]] && q->N[a] < q->N[b]);
}

static int compare_double(const struct tsqubo_queue *q, const void *v_, size_t a, size_t b) {
  const double *v = (const double *)v_;
  return v[q->N[a]] < v[q->N[b]] || (v[q->N[a]] == v[q->N[b]] && q->N[a] < q->N[b]);
}

static void queue_heapify(struct tsqubo_queue *q, const void *dx,
                          int (*compare)(const struct tsqubo_queue *, const void *, size_t, size_t),
                          size_t k) {
  for (;;) {
    size_t smallest = k, left = 2 * k + 1, right = 2 * k + 2;
    if (left < q->size && compare(q, dx, left, k)) smallest = left;
    if (right < q->size && compare(q, dx, right, smallest)) smallest = right;
    if (smallest == k) break;
    queue_reorder(q, q->N[k], smallest);
    k = smallest;
  }
}

static void queue_decrease(struct tsqubo_queue *q, const void *dx,
                           int (*compare)(const struct tsqubo_queue *, const void *, size_t,
                                          size_t),
                           size_t k) {
  for (size_t parent = (k - 1) / 2; k && compare(q, dx, k, parent); parent = (k - 1) / 2) {
    queue_reorder(q, q->N[k], parent);
    k = parent;
  }
}

static void queue_push(struct tsqubo_queue *q, const void *dx,
                       int (*compare)(const struct tsqubo_queue *, const void *, size_t, size_t),
                       size_t i) {
  size_t k = q->size++;
  queue_reorder(q, i, k);
  queue_decrease(q, dx, compare, k);
}

static void queue_pop(struct tsqubo_queue *q, const void *dx,
                      int (*compare)(const struct tsqubo_queue *, const void *, size_t, size_t)) {
  queue_reorder(q, queue_top(q), --q->size);
  queue_heapify(q, dx, compare, 0);
}

static void queue_remove(struct tsqubo_queue *q, const void *dx,
                         int (*compare)(const struct tsqubo_queue *, const void *, size_t, size_t),
                         size_t k) {
  for (size_t parent = (k - 1) / 2; k; parent = (k - 1) / 2) {
    queue_reorder(q, q->N[k], parent);
    k = parent;
  }
  queue_pop(q, dx, compare);
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
  size_t *tabulist;
  /** The iteration counter. */
  size_t iteration;
  /** A priority queue of indices not in the tabu list, based on reevaluation vector values. */
  struct tsqubo_queue d;
  /** A priority queue of indices in the tabu list, based on the values in the tabu list. */
  struct tsqubo_queue l;
  /** A priority queue of indices in the tabu list, based on reevaluation vector values. */
  struct tsqubo_queue ld;
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
  ts->tabulist = (size_t *)calloc(inst->n, sizeof(size_t));
  queue_init(&ts->d, inst->n);
  queue_init(&ts->l, inst->n);
  queue_init(&ts->ld, inst->n);
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
  free(ts->tabulist);
  queue_free(&ts->d);
  queue_free(&ts->l);
  queue_free(&ts->ld);
}

/**
 * Copies the current solution to the incumbent solution.
 *
 * @param ts an initialized TS structure.
 */
void tsqubo_commit_incumbent(struct tsqubo *ts) {
  ts->inc.fx = ts->cur.fx;
  memcpy(ts->inc.x, ts->cur.x, ts->inst.n * sizeof(double));
  memcpy(ts->inc.dx, ts->cur.dx, ts->inst.n * sizeof(double));
}

/**
 * Zeroes the current and incumbent solution.
 *
 * @param ts an initialized TS structure.
 */
void tsqubo_reset_solutions(struct tsqubo *ts) {
  ts->cur.fx = 0;
  memset(ts->cur.x, 0, ts->inst.n * sizeof(double));
  for (size_t i = 0; i < ts->inst.n; i++) {
    ts->cur.dx[i] = 0;
    for (size_t k = ts->inst.R[i]; k < ts->inst.R[i + 1]; k++)
      if (ts->inst.C[k] == i) {
        ts->cur.dx[i] = ts->inst.Q[k];
        break;
      }
  }
  tsqubo_commit_incumbent(ts);
}

/**
 * Resets the tabu list.
 *
 * @param ts an initialized TS structure.
 */
void tsqubo_reset_tabu(struct tsqubo *ts) {
  memset(ts->tabulist, 0, ts->inst.n * sizeof(size_t));
  ts->iteration = 0;
  ts->d.size = ts->l.size = ts->ld.size = 0;
  for (size_t i = 0; i < ts->inst.n; i++) queue_push(&ts->d, ts->cur.dx, compare_double, i);
}

/**
 * Flips a variable of the current solution.
 *
 * @param ts an initialized TS structure.
 * @param i the index of the variable to flip.
 */
void tsqubo_flip_current(struct tsqubo *ts, size_t i) {
  ts->cur.fx += ts->cur.dx[i];
  ts->cur.x[i] = 1 - ts->cur.x[i];
  for (size_t k = ts->inst.R[i]; k < ts->inst.R[i + 1]; k++) {
    size_t j = ts->inst.C[k];
    double d = ts->cur.dx[j];
    ts->cur.dx[j] =
        j == i ? -ts->cur.dx[j]
               : ts->cur.dx[j] - (1 - 2 * ts->cur.x[i]) * (1 - 2 * ts->cur.x[j]) * ts->inst.Q[k];
    if (ts->cur.dx[j] < d) {
      if (queue_contains(&ts->d, j)) queue_decrease(&ts->d, ts->cur.dx, compare_double, ts->d.I[j]);
      if (queue_contains(&ts->ld, j))
        queue_decrease(&ts->ld, ts->cur.dx, compare_double, ts->ld.I[j]);
    }
    if (ts->cur.dx[j] > d) {
      if (queue_contains(&ts->d, j)) queue_heapify(&ts->d, ts->cur.dx, compare_double, ts->d.I[j]);
      if (queue_contains(&ts->ld, j))
        queue_heapify(&ts->ld, ts->cur.dx, compare_double, ts->ld.I[j]);
    }
  }
}

/**
 * Performs a local search on the incumbent solution of the search.
 *
 * @param ts an initialized TS structure.
 */
void tsqubo_local_search(struct tsqubo *ts) {
  for (;;) {
    size_t i = queue_top(&ts->d);
    if (ts->ld.size) {
      size_t j = queue_top(&ts->ld);
      if (ts->cur.dx[j] < ts->cur.dx[i]) i = j;
    }
    if (ts->cur.dx[i] >= 0) break;
    tsqubo_flip_current(ts, i);
  }
}

/**
 * Performs a TS iteration.
 *
 * @param ts an initialized TS structure.
 * @param K the tabu tenure to assign to flipped variables.
 * @return whether an improved solution was found or not.
 */
int tsqubo_iterate(struct tsqubo *ts, size_t K) {
  size_t i = queue_top(&ts->d);
  if (ts->ld.size) {
    size_t j = queue_top(&ts->ld);
    if (ts->cur.fx + ts->cur.dx[j] < ts->inc.fx && ts->cur.dx[j] < ts->cur.dx[i]) i = j;
  }
  ++ts->iteration;
  ts->tabulist[i] = ts->iteration + K;
  if (queue_contains(&ts->l, i)) {
    queue_heapify(&ts->l, ts->tabulist, compare_size, ts->l.I[i]);
  } else {
    queue_pop(&ts->d, ts->cur.dx, compare_double);
    queue_push(&ts->l, ts->tabulist, compare_size, i);
    queue_push(&ts->ld, ts->cur.dx, compare_double, i);
  }
  while (ts->l.size && ts->iteration >= ts->tabulist[queue_top(&ts->l)]) {
    size_t j = queue_top(&ts->l);
    queue_pop(&ts->l, ts->tabulist, compare_size);
    queue_remove(&ts->ld, ts->cur.dx, compare_double, ts->ld.I[j]);
    queue_push(&ts->d, ts->cur.dx, compare_double, j);
  }
  tsqubo_flip_current(ts, i);
  return ts->cur.fx < ts->inc.fx;
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
  for (size_t i = 0; i < cutoff; i++)
    if (tsqubo_iterate(ts, K)) {
      tsqubo_local_search(ts);
      tsqubo_commit_incumbent(ts);
      goto begin;
    }
}