# libtsqubo

A header-only C library for solving QUBO problems through Tabu Search.

## Basic usage

```c
#include "tsqubo.h"

int main() {
  struct tsqubo_instance* instance = tsqubo_instance_new();
  tsqubo_instance_add_component(instance, 0, 0, 2);
  tsqubo_instance_add_component(instance, 1, 1, -3);
  tsqubo_instance_add_component(instance, 2, 2, 5);
  tsqubo_instance_add_component(instance, 0, 1, -1);

  struct tsqubo* ts = tsqubo_new(instance);
  tsqubo_instance_free(instance);
  free(instance);

  size_t tabu_tenure_constant = 10, improvement_cutoff = 5;
  tsqubo_iterate_cutoff(ts, tabu_tenure_constant, improvement_cutoff);
  printf("%lf\n", ts->inc.fx);

  tsqubo_free(ts);
  free(ts);
  return 0;
}
```
