#include <assert.h>
#include <stdio.h>

#include "tsqubo.h"

int main() {
  FILE* file = fopen("G1", "r");
  size_t n, nonzeros;
  assert(fscanf(file, "%zu%zu", &n, &nonzeros) == 2);
  assert(n > 0);
  assert(nonzeros > 0);
  double* Q = calloc(n * n, sizeof(double));
  struct tsqubo_instance inst;
  tsqubo_instance_init(&inst, 3);
  double* diag = calloc(n, sizeof(double));
  for (size_t k = 0; k < nonzeros; k++) {
    size_t i, j;
    double q;
    assert(fscanf(file, "%zu%zu%lf", &i, &j, &q) == 3);
    i--;
    j--;
    Q[n * i + i] -= q;
    Q[n * j + j] -= q;
    Q[n * i + j] += 2 * q;
    Q[n * j + i] += 2 * q;
    diag[i] -= q;
    diag[j] -= q;
    tsqubo_instance_add_component(&inst, i, j, q);
  }
  fclose(file);
  for (size_t i = 0; i < n; i++) tsqubo_instance_add_component(&inst, i, i, diag[i]);
  free(diag);
  struct tsqubo* ts = tsqubo_new(&inst);
  tsqubo_instance_free(&inst);
  for (size_t i = 0; i < 10; i++) {
    tsqubo_reset_solutions(ts);
    tsqubo_reset_tabu(ts);
    for (size_t j = 0; j < n; j++)
      if (rand() % 2) {
        tsqubo_flip_current(ts, j);
        double fy = 0;
        for (size_t i = 0; i < n; i++)
          for (size_t j = 0; j <= i; j++) fy += ts->cur.x[i] * ts->cur.x[j] * Q[n * i + j];
        assert(fy == ts->cur.fx);
      }
    tsqubo_commit_incumbent(ts);
    double fx = 0;
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j <= i; j++) fx += ts->inc.x[i] * ts->inc.x[j] * Q[n * i + j];
    assert(fx == ts->inc.fx);
    for (size_t l = 0; l < 100; l++) {
      if (tsqubo_iterate(ts, n / 10)) {
        tsqubo_local_search(ts);
        tsqubo_commit_incumbent(ts);
      }
      double fx = 0;
      for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j <= i; j++) fx += ts->inc.x[i] * ts->inc.x[j] * Q[n * i + j];
      assert(fx == ts->inc.fx);
      double fy = 0;
      for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j <= i; j++) fy += ts->cur.x[i] * ts->cur.x[j] * Q[n * i + j];
      assert(fy == ts->cur.fx);
    }
  }
  tsqubo_free(ts);
  free(ts);
  free(Q);
  return 0;
}