#pragma once

#include <stdbool.h>
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
 * Initializes a solution structure.
 *
 * @param sol pointer to an uninitialized solution structure.
 * @param n the number of variables.
 */
void tsqubo_solution_init(struct tsqubo_solution *sol, size_t n) {
  sol->fx = 0;
  sol->x = (double *)calloc(n, sizeof(double));
  sol->dx = (double *)calloc(n, sizeof(double));
}

/**
 * Frees a solution structure.
 *
 * @param sol pointer to an initialized solution structure.
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

/** A QUBO problem instance structure in CSR format. */
struct tsqubo_csr_instance {
  /** The number of variables, or rows in @ref Q. */
  size_t n;
  /** The coefficients of the matrix. */
  double *Q;
  /** The columns of rows. */
  size_t *C;
  /** The position in @ref cols containing the columns of each row. */
  size_t *R;
};

static void tsqubo_csr_instance_init(struct tsqubo_csr_instance *cinst,
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
}

static void tsqubo_csr_instance_free(struct tsqubo_csr_instance *cinst) {
  free(cinst->Q);
  free(cinst->C);
  free(cinst->R);
}

#ifdef TSQUBO_SPARSE

/** An indexed priority queue structure. */
struct ipq {
  /** A permutation of `{1,2,...n}`. */
  size_t *N;
  /** An index vector that associates each `{1,2,...n}` to its position in @ref N. */
  size_t *I;
  /** The size of the queue. */
  size_t size;
};

static void ipq_init(struct ipq *q, size_t n) {
  q->N = (size_t *)calloc(n, sizeof(size_t));
  q->I = (size_t *)calloc(n, sizeof(size_t));
  for (size_t i = 0; i < n; i++) q->N[i] = q->I[i] = i;
  q->size = 0;
}

static void ipq_free(struct ipq *q) {
  free(q->N);
  free(q->I);
}

static void ipq_reorder(struct ipq *q, size_t i, size_t k) {
  q->I[q->N[k]] = q->I[i];
  q->N[q->I[i]] = q->N[k];
  q->I[i] = k;
  q->N[k] = i;
}

static int ipq_contains(const struct ipq *q, size_t i) { return q->I[i] < q->size; }

static size_t ipq_top(const struct ipq *q) { return q->N[0]; }

static int compare_size(const struct ipq *q, const void *v_, size_t a, size_t b) {
  const size_t *v = (const size_t *)v_;
  return v[q->N[a]] < v[q->N[b]] || (v[q->N[a]] == v[q->N[b]] && q->N[a] < q->N[b]);
}

static int compare_double(const struct ipq *q, const void *v_, size_t a, size_t b) {
  const double *v = (const double *)v_;
  return v[q->N[a]] < v[q->N[b]] || (v[q->N[a]] == v[q->N[b]] && q->N[a] < q->N[b]);
}

static void ipq_sift_up(struct ipq *q, const void *dx,
                        int (*compare)(const struct ipq *, const void *, size_t, size_t),
                        size_t k) {
  for (;;) {
    size_t smallest = k, left = 2 * k + 1, right = 2 * k + 2;
    if (left < q->size && compare(q, dx, left, k)) smallest = left;
    if (right < q->size && compare(q, dx, right, smallest)) smallest = right;
    if (smallest == k) break;
    ipq_reorder(q, q->N[k], smallest);
    k = smallest;
  }
}

static void ipq_sift_down(struct ipq *q, const void *dx,
                          int (*compare)(const struct ipq *, const void *, size_t, size_t),
                          size_t k) {
  for (size_t parent = (k - 1) / 2; k && compare(q, dx, k, parent); parent = (k - 1) / 2) {
    ipq_reorder(q, q->N[k], parent);
    k = parent;
  }
}

static void ipq_push(struct ipq *q, const void *dx,
                     int (*compare)(const struct ipq *, const void *, size_t, size_t), size_t i) {
  size_t k = q->size++;
  ipq_reorder(q, i, k);
  ipq_sift_down(q, dx, compare, k);
}

static void ipq_pop(struct ipq *q, const void *dx,
                    int (*compare)(const struct ipq *, const void *, size_t, size_t)) {
  ipq_reorder(q, ipq_top(q), --q->size);
  ipq_sift_up(q, dx, compare, 0);
}

static void ipq_remove(struct ipq *q, const void *dx,
                       int (*compare)(const struct ipq *, const void *, size_t, size_t), size_t k) {
  for (size_t parent = (k - 1) / 2; k; parent = (k - 1) / 2) {
    ipq_reorder(q, q->N[k], parent);
    k = parent;
  }
  ipq_pop(q, dx, compare);
}

#endif

/** A QUBO tabu search structure. */
struct tsqubo {
  /** The problem instance. */
  struct tsqubo_csr_instance inst;
  /** The incumbent solution of the search. */
  struct tsqubo_solution inc;
  /** The current solution of the search. */
  struct tsqubo_solution cur;
  /** The tabu list. */
  size_t *tabulist;
  /** The iteration counter. */
  size_t iteration;
#ifdef TSQUBO_SPARSE
  /** A priority queue of indices not in the tabu list, based on reevaluation vector values. */
  struct ipq d;
  /** A priority queue of indices in the tabu list, based on the values in the tabu list. */
  struct ipq l;
  /** A priority queue of indices in the tabu list, based on reevaluation vector values. */
  struct ipq ld;
#endif
};

/**
 * Initializes the TS structure.
 *
 * @param ts an uninitialized TS structure.
 */
void tsqubo_init(struct tsqubo *ts, struct tsqubo_instance *inst) {
  tsqubo_csr_instance_init(&ts->inst, inst);
  tsqubo_solution_init(&ts->inc, inst->n);
  tsqubo_solution_init(&ts->cur, inst->n);
  ts->tabulist = (size_t *)calloc(inst->n, sizeof(size_t));
#ifdef TSQUBO_SPARSE
  ipq_init(&ts->d, inst->n);
  ipq_init(&ts->l, inst->n);
  ipq_init(&ts->ld, inst->n);
#endif
}

/**
 * Allocates and initializes the TS structure.
 *
 * @param inst the problem instance to solve.
 * @return a newly-allocated TS structure.
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
  tsqubo_csr_instance_free(&ts->inst);
  tsqubo_solution_free(&ts->inc);
  tsqubo_solution_free(&ts->cur);
  free(ts->tabulist);
#ifdef TSQUBO_SPARSE
  ipq_free(&ts->d);
  ipq_free(&ts->l);
  ipq_free(&ts->ld);
#endif
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
#ifdef TSQUBO_SPARSE
  ts->d.size = ts->l.size = ts->ld.size = 0;
  for (size_t i = 0; i < ts->inst.n; i++) ipq_push(&ts->d, ts->cur.dx, compare_double, i);
#endif
}

/**
 * Flips a variable of the current solution to its complementary value.
 *
 * @param ts an initialized TS structure.
 * @param i the index of the variable to flip.
 */
void tsqubo_flip_current(struct tsqubo *ts, size_t i) {
  ts->cur.fx += ts->cur.dx[i];
  ts->cur.x[i] = 1 - ts->cur.x[i];
  ts->cur.dx[i] *= -1;
  for (size_t k = ts->inst.R[i]; k < ts->inst.R[i + 1]; k++) {
    size_t j = ts->inst.C[k];
#ifdef TSQUBO_SPARSE
    double d = ts->cur.dx[j];
#endif
    if (i != j)
    ts->cur.dx[j] -= (1 - 2 * ts->cur.x[i]) * (1 - 2 * ts->cur.x[j]) * ts->inst.Q[k];
#ifdef TSQUBO_SPARSE
    if (ts->cur.dx[j] < d) {
      if (ipq_contains(&ts->d, j)) ipq_sift_down(&ts->d, ts->cur.dx, compare_double, ts->d.I[j]);
      if (ipq_contains(&ts->ld, j)) ipq_sift_down(&ts->ld, ts->cur.dx, compare_double, ts->ld.I[j]);
    }
    if (ts->cur.dx[j] > d) {
      if (ipq_contains(&ts->d, j)) ipq_sift_up(&ts->d, ts->cur.dx, compare_double, ts->d.I[j]);
      if (ipq_contains(&ts->ld, j)) ipq_sift_up(&ts->ld, ts->cur.dx, compare_double, ts->ld.I[j]);
    }
#endif
  }
}

/**
 * Performs a local search on the incumbent solution of the search.
 *
 * @param ts an initialized TS structure.
 */
void tsqubo_local_search(struct tsqubo *ts) {
  for (;;) {
#ifdef TSQUBO_SPARSE
    size_t i = ipq_top(&ts->d);
    if (ts->ld.size) {
      size_t j = ipq_top(&ts->ld);
      if (ts->cur.dx[j] < ts->cur.dx[i]) i = j;
    }
#else
    size_t i = 0;
    for (size_t j = 1; j < ts->inst.n; j++)
      if (ts->cur.dx[j] < ts->cur.dx[i]) i = j;
#endif
    if (ts->cur.dx[i] >= 0) break;
    tsqubo_flip_current(ts, i);
  }
}

/**
 * Performs a TS iteration.
 *
 * @param ts an initialized TS structure.
 * @param ttc the tabu tenure to assign to flipped variables.
 * @return whether an improved solution was found or not.
 */
bool tsqubo_iterate(struct tsqubo *ts, size_t ttc) {
  size_t i = ts->inst.n;
#ifdef TSQUBO_SPARSE
  if (ts->d.size > 0) {
    i = ipq_top(&ts->d);
  }
  if (ts->ld.size) {
    size_t j = ipq_top(&ts->ld);
    if (ts->cur.fx + ts->cur.dx[j] < ts->inc.fx && ts->cur.dx[j] < ts->cur.dx[i]) i = j;
  }
#else
  for (size_t j = 0; j < ts->inst.n; j++) {
    if (ts->iteration < ts->tabulist[j] && ts->cur.fx + ts->cur.dx[j] >= ts->inc.fx) continue;
    if (i == ts->inst.n || ts->cur.dx[j] < ts->cur.dx[i]) i = j;
  }
#endif
  ++ts->iteration;
  if (i == ts->inst.n) 
  // If all flips are tabu and none improve the incumbent, don't flip.
  return 0;

  ts->tabulist[i] = ts->iteration + ttc;
#ifdef TSQUBO_SPARSE
  if (ipq_contains(&ts->l, i)) {
    ipq_sift_up(&ts->l, ts->tabulist, compare_size, ts->l.I[i]);
  } else {
    ipq_pop(&ts->d, ts->cur.dx, compare_double);
    ipq_push(&ts->l, ts->tabulist, compare_size, i);
    ipq_push(&ts->ld, ts->cur.dx, compare_double, i);
  }
  while (ts->l.size && ts->iteration >= ts->tabulist[ipq_top(&ts->l)]) {
    size_t j = ipq_top(&ts->l);
    ipq_pop(&ts->l, ts->tabulist, compare_size);
    ipq_remove(&ts->ld, ts->cur.dx, compare_double, ts->ld.I[j]);
    ipq_push(&ts->d, ts->cur.dx, compare_double, j);
  }
#endif
  tsqubo_flip_current(ts, i);
  return ts->cur.fx < ts->inc.fx;
}

/**
 * Performs TS iterations until no improvements are possible for a number of iterations.
 *
 * @param ts an initialized TS structure.
 * @param ttc the tabu tenure to assign to flipped variables.
 * @param cutoff the improvement cutoff.
 */
void tsqubo_iterate_cutoff(struct tsqubo *ts, size_t ttc, size_t cutoff) {
begin:
  for (size_t i = 0; i < cutoff; i++)
    if (tsqubo_iterate(ts, ttc)) {
      tsqubo_local_search(ts);
      tsqubo_commit_incumbent(ts);
      goto begin;
    }
}