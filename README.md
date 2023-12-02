# Genetic Algorithm genetic algorithm for the binary knapsack problem.

C++ implementation of a genetic algorithm for the binary knapsack problem.

### How to use:
- Compile with:
```bash
g++ -O3 genetic_algo_knapsack.cpp -o genetic_algo_knapsack
```

- Usage:
Flags:
- -r : Binary number to enable or disable using random items for the knapsack.
- -g : Integer for the number of generations the algorithm is going to run.
- -s : Integer for the number of solutions per population for the algorithm.

```bash
./genetic_algo_knapsack -r 1 -g 10000 -s 20
```

(**Note:** If random is disabled, you must provide the itens for the knapsack in the format (weight value) via stdin).
