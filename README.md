# libtsqubo

A header-only C/C++ library for solving QUBO problems through Tabu Search.

# Features

* The instance matrix is stored in CSR format;
* Using the compile-time flag `TSQUBO_SPARSE` makes each Tabu Search iteration `O(D log(n))`-time,
where `D` is the maximum variable degree.

# Installation

Simply download the `tsqubo.h` file and include it in your source.

To run the sanity tests, clone the repository and run the command `make`.

# Example usage

```c
#define TSQUBO_SPARSE
#include "tsqubo.h"

int main() {
  // create a QUBO instance and add some components to the matrix
  struct tsqubo_instance *instance = tsqubo_instance_new();
  tsqubo_instance_add_component(instance, 0, 0, 2);
  tsqubo_instance_add_component(instance, 1, 1, -3);
  tsqubo_instance_add_component(instance, 2, 2, 5);
  tsqubo_instance_add_component(instance, 0, 1, -1);

  // create a TS instance that will solve the QUBO instance
  struct tsqubo *ts = tsqubo_new(instance);
  tsqubo_instance_free(instance);
  free(instance);

  // run TS until no improved solutions are found for 5 consecutive iterations
  size_t tabu_tenure_constant = 10, improvement_cutoff = 5;
  tsqubo_iterate_cutoff(ts, tabu_tenure_constant, improvement_cutoff);

  // print the best solution value
  printf("%lf\n", ts->inc.fx);

  tsqubo_free(ts);
  free(ts);
  return 0;
}
```

# API

Here we describe the structures and functions that users may use to solve QUBO problem instances by Tabu Search.

## Instance manipulation

The following structures and functions can be used to initialize and manipulate a QUBO problem instance.

### `struct tsqubo_instance`
Represents a QUBO problem instance structure.
Users can use the functions below to initialize a problem instance to be solved by TS.

### `struct tsqubo_instance *tsqubo_instance_new(size_t capacity)`
Allocates and initializes a problem instance structure.

Parameters:
* `capacity`: an initial capacity for the number of nonzero components in the matrix.

Returns:
A newly-allocated QUBO problem instance structure.

### `void tsqubo_instance_free(struct tsqubo_instance *inst)`
Frees the memory associated to problem instance structure.

Parameters:
* `inst`: an initialized problem instance structure.

### `void tsqubo_instance_add_component(struct tsqubo_instance *inst, size_t i, size_t j, double q)`
Sets a coefficient of the instance matrix.

Parameters:
* `inst`: an initialized problem instance structure.
* `i`: the row of the matrix.
* `j`: the column of the matrix.
* `q`: the value of the coefficient.

## Tabu Search conduction

The following structures and functions contain information pertaining to a Tabu Search procedure, and can be used to conduct the search.

### `struct tsqubo_solution`
A QUBO solution structure.

Fields:
* `double fx`: The objective function value.
* `double *x`: The solution vector.

### `struct tsqubo`
A QUBO tabu search structure.

Fields:
* `struct tsqubo_solution inc`: The incumbent solution of the search.
* `struct tsqubo_solution cur`: The current solution of the search.

### `struct tsqubo *tsqubo_new(struct tsqubo_instance *inst)`
Allocates and initializes a TS structure.

Parameters:
* `inst`: the problem instance to solve.

Returns:
A newly-allocated TS structure.

### `void tsqubo_free(struct tsqubo *ts)`
Frees memory associated with a TS structure.

Parameters:
* `ts`: an initialized TS structure.

### `void tsqubo_reset_tabu(struct tsqubo *ts)`
Resets the tabu list.

Parameters:
* `ts`: an initialized TS structure.

### `void tsqubo_reset_solutions(struct tsqubo *ts)`
Resets the current and incumbent solutions, setting their components to zeroes.

Parameters:
* `ts`: an initialized TS structure.

### `void tsqubo_flip_current(struct tsqubo *ts, size_t i)`
Flips a variable of the current solution to its complementary value.

Parameters:
* `ts`: an initialized TS structure.
* `i`: the index of the variable to flip.

### `void tsqubo_commit_incumbent(struct tsqubo *ts)`
Copies the current solution to the incumbent solution.

Parameters:
* `ts`: an initialized TS structure.

### `bool tsqubo_iterate(struct tsqubo *ts, size_t ttc)`
Performs a TS iteration.

Parameters:
* `ts`: an initialized TS structure.
* `ttc`: the tabu tenure to assign to flipped variables.

Returns:
Whether an improved solution was found or not.

### `void tsqubo_iterate_cutoff(struct tsqubo *ts, size_t ttc, size_t cutoff)`
Performs TS iterations until no improvements are possible for a number of iterations.

Parameters:
* `ts`: an initialized TS structure.
* `ttc`: the tabu tenure to assign to flipped variables.
* `cutoff`: the improvement cutoff.

## Internal structures and functions

These structures and functions are used internally and should not need to be used by users explicitly.

### `void tsqubo_instance_init(struct tsqubo_instance *inst, size_t capacity)`
Initializes the problem instance structure.
 
Parameters:
* `inst`: an uninitialized problem instance structure.
* `capacity`: the initial capacity for components in the matrix.

### `void tsqubo_solution_init(struct tsqubo_solution *sol, size_t n)`
Initializes a solution structure.
 
Parameters:
* `sol`: an uninitialized solution structure.
* `n`: the number of variables.

### `void tsqubo_solution_free(struct tsqubo_solution *sol)`
Frees the memory associated with a solution structure.

Parameters:
* `sol`: pointer to an initialized solution structure.

### `struct tsqubo_csr_instance`
A QUBO problem instance structure in CSR format.

### `void tsqubo_csr_instance_init(struct tsqubo_csr_instance *cinst, struct tsqubo_instance *inst)`
Converts a QUBO problem instance to CSR format.

Parameters:
* `cinst` pointer to an uninitialized QUBO problem instance in CSR format.
* `inst` pointer to an initialized QUBO problem instance.

### `void tsqubo_csr_instance_free(struct tsqubo_csr_instance *cinst)`
Frees the memory associated with a QUBO problem instance in CSR format.

Parameters:
* `cinst` pointer to an initialized QUBO problem instance in CSR format.

### `struct queue`
An indexed priority queue structure. Enabled with the compile-time flag `TSQUBO_SPARSE`.

### `void queue_init(struct queue *q, size_t n)`
Initializes an indexed priority queue data structure.

Parameters:
* `q`: an uninitialized indexed priority queue data structure.
* `n`: states that the structure may store values from the set `{1,2,...n}`.

### `void queue_free(struct queue *q)`
Frees the memory associated with an indexed priority queue structure.

Parameters:
* `q`: an initialized indexed priority queue data structure.

### `void queue_reorder(struct queue *q, size_t i, size_t k)`
Changes the position of an element within the queue.

Parameters:
* `q`: an initialized indexed priority queue data structure.
* `i`: The element from `{1,2,...n}` which to change the position of.
* `k`: The new position of `i` in the queue.

### `size_t queue_top(const struct queue *q)`
Obtains the element from `{1,2,...,n}` with the minimum priority in the queue.

Parameters:
* `q`: an initialized indexed priority queue data structure.

Returns:
An element from `{1,2,...,n}`.

### `void queue_heapify(struct queue *q, const void *dx, int (*compare)(const struct queue *, const void *, size_t, size_t), size_t k)`
Corrects the minimum-heap property starting at a given position up to the end of the queue.

Parameters:
* `q`: an initialized indexed priority queue data structure.
* `dx`: an array where each `i`th element is the priority of the element `i`.
* `compare`: an appropriate comparison function according to the type of `dx`.
* `k`: The starting index from which to correct the minimum-heap property.

### `void queue_decrease(struct queue *q, const void *dx, int (*compare)(const struct queue *, const void *, size_t, size_t), size_t k)`
Corrects the minimum-heap property starting at a given position down to the beginning of the queue.

Parameters:
* `q`: an initialized indexed priority queue data structure.
* `dx`: an array where each `i`th element is the priority of the element `i`.
* `compare`: an appropriate comparison function according to the type of `dx`.
* `k`: The end index up to which to correct the minimum-heap property.

### `void queue_push(struct queue *q, const void *dx, int (*compare)(const struct queue *, const void *, size_t, size_t), size_t i)`
Inserts an element from `{1,2,...n}` into the queue.

Parameters:
* `q`: an initialized indexed priority queue data structure.
* `dx`: an array where each `i`th element is the priority of the element `i`.
* `compare`: an appropriate comparison function according to the type of `dx`.
* `i`: the element to insert into the queue.

### `void queue_pop(struct queue *q, const void *dx, int (*compare)(const struct queue *, const void *, size_t, size_t))`
Removes the top element of the queue.

Parameters:
* `q`: an initialized indexed priority queue data structure.
* `dx`: an array where each `i`th element is the priority of the element `i`.
* `compare`: an appropriate comparison function according to the type of `dx`.

### `void queue_remove(struct queue *q, const void *dx, int (*compare)(const struct queue *, const void *, size_t, size_t), size_t k)`
Removes an element at an arbitrary position in the queue.

Parameters:
* `q`: an initialized indexed priority queue data structure.
* `dx`: an array where each `i`th element is the priority of the element `i`.
* `compare`: an appropriate comparison function according to the type of `dx`.
* `k`: the position of the element to remove.

### `int compare_size(const struct tsqubo_queue *q, const void *v_, size_t a, size_t b)`
Comparison function used in `queue_heapify` and `queue_decrease`, when the key values are of type `size_t`.

### `int compare_double(const struct tsqubo_queue *q, const void *v_, size_t a, size_t b)`
Comparison function used in `queue_heapify` and `queue_decrease`, when the key values are of type `double`.

### `int compare_instance_components(const void *a_, const void *b_)`
Comparison function used to sort the components in a QUBO problem instance before converting it to CSR format.