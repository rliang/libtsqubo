# libtsqubo

A header-only C/C++ library for solving QUBO problems through Tabu Search.

## Features

* The instance matrix is stored in CSR format;
* Using the compile-time flag `TSQUBO_SPARSE` makes each Tabu Search iteration `O(D log(n))`-time,
where `D` is the maximum variable degree.

## Example usage

```c
#include "tsqubo.h"

#define TSQUBO_SPARSE

int main() {
  // create a QUBO instance and add some components to the matrix
  struct tsqubo_instance* instance = tsqubo_instance_new();
  tsqubo_instance_add_component(instance, 0, 0, 2);
  tsqubo_instance_add_component(instance, 1, 1, -3);
  tsqubo_instance_add_component(instance, 2, 2, 5);
  tsqubo_instance_add_component(instance, 0, 1, -1);

  // create a TS instance that will solve the QUBO instance
  struct tsqubo* ts = tsqubo_new(instance);
  tsqubo_instance_free(instance);
  free(instance);

  // run TS until no improved solutions are found for `5 * n` consecutive iterations
  size_t tabu_tenure_constant = 10, improvement_cutoff = 5;
  tsqubo_iterate_cutoff(ts, tabu_tenure_constant, improvement_cutoff);

  // print the best solution value
  printf("%lf\n", ts->inc.fx);

  tsqubo_free(ts);
  free(ts);
  return 0;
}
```

## API

## Internal API

## Running tests

Simply run the command:
```sh
make
```
