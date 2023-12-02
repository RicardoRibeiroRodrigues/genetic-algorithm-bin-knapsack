#include <random>
#include <vector>
#include <iostream>
#include <chrono>
#include <string>

using namespace std;

struct Item
{
    int id;
    int weight;
    int value;
};
typedef std::chrono::high_resolution_clock myclock;
myclock::time_point beginning = myclock::now();

typedef vector<vector<int>> population;

void print_population(population pop)
{
    cout << '[' << endl;
    for (int i = 0; i < pop.size(); i++)
    {
        cout << '[';
        for (auto value : pop[i])
            cout << value << ' ';
        cout << ']' << endl;
    }
    cout << ']' << endl;
}

void print_vec(vector<int> vec)
{
    for (auto item : vec)
        cout << item << ' ';
    cout << endl;
}

// Calculate the fitness value of each solution in the current population.
vector<int> fitness(vector<int> weight, vector<int> value, population pop, int capacity)
{
    vector<int> fit(pop.size());

    for (int i = 0; i < pop.size(); i++)
    {
        int sum_w = 0, sum_v = 0;
        vector<int> curr_ind = pop[i];
        for (int j = 0; j < curr_ind.size(); j++)
        {
            sum_w += curr_ind[j] * weight[j];
            sum_v += curr_ind[j] * value[j];
        }

        if (sum_w <= capacity)
            fit[i] = sum_v;
        else
            // fit[i] = 0;
            // Resultado melhor:
            fit[i] = capacity - sum_w;
    }

    return fit;
}

// Return the index of the maximum value in the numeric array.
int idx_max(vector<int> numeric_array)
{
    int max_value = numeric_array[0];
    int idx_max = 0;
    for (int i = 1; i < numeric_array.size(); i++)
    {
        if (numeric_array[i] > max_value)
        {
            max_value = numeric_array[i];
            idx_max = i;
        }
    }
    return idx_max;
}

// Select the best individuals in the current generation as parents for producing the offspring of the next generation.
population selection(vector<int> fitness, int num_parents, population pop)
{
    population parents;
    parents.reserve(num_parents);

    for (int i = 0; i < num_parents; i++)
    {
        int max_fit_idx = idx_max(fitness);
        parents.push_back(pop[max_fit_idx]);
        // Set the fitness of the current best to -9999 so that it is not selected as parent again.
        fitness[max_fit_idx] = -9999;
    }
    return parents;
}

// Implement crossover in the offspring.
population crossover(population parents, int num_offspring)
{
    population offsprings;
    offsprings.reserve(num_offspring);

    int crossover_point = parents[0].size() / 2;
    const float crossover_rate = 0.8;
    int i = 0;
    int cnt_offspring = 0;

    // Random double from 0 to 1.
    myclock::duration d = myclock::now() - beginning;
    unsigned seed = d.count();
    default_random_engine generator(seed);
    uniform_real_distribution<double> distribution(0.0, 1.0);

    while (cnt_offspring < num_offspring)
    {
        int parent_1_idx = i % parents.size();
        int parent_2_idx = (i + 1) % parents.size();

        double x = distribution(generator);
        if (x <= crossover_rate)
        {
            vector<int> parent1 = parents[parent_1_idx];
            vector<int> parent2 = parents[parent_2_idx];

            vector<int> cross_res(parent1.size());
            // Copy from 0...crossover_point to res.
            copy(parent1.begin(), parent1.begin() + crossover_point - 1, cross_res.begin());
            // Copy from crossover_point...end to res.
            copy(parent2.begin() + crossover_point, parent2.end(), cross_res.begin() + crossover_point);

            offsprings.push_back(cross_res);
            cnt_offspring += 1;
        }
        i += 1;
    }
    return offsprings;
}

// Implement mutation in the offspring.
population mutation(population offspring)
{
    population mutants = offspring;
    const float mutation_rate = 0.15;

    // Random double from 0 to 1.
    myclock::duration d = myclock::now() - beginning;
    unsigned seed = d.count();
    default_random_engine generator(seed);
    uniform_real_distribution<double> distribution(0.0, 1.0);
    uniform_int_distribution<int> distri_int(0, offspring[0].size() - 1);

    for (int i = 0; i < offspring.size(); i++)
    {
        double random_value = distribution(generator);

        if (random_value > mutation_rate)
            continue;

        int int_random_value = distri_int(generator);
        // Flip the bit;
        mutants[i][int_random_value] = !mutants[i][int_random_value];
    }
    return mutants;
}

// Execute the genetic algorithm.
vector<int> optimize(vector<int> weight, vector<int> value, population &pop, int pop_size, int num_generations, int capacity)
{
    int num_parents = pop_size / 2;
    int num_offspring = pop_size - num_parents;

    for (int i = 0; i < num_generations; i++)
    {
        vector<int> fit = fitness(weight, value, pop, capacity);
        population parents = selection(fit, num_parents, pop);
        population offsprings = crossover(parents, num_offspring);
        population mutants = mutation(offsprings);

        // cout << "FIT " << i << endl;
        // print_vec(fit);

        // Copy from 0...crossover_point to res.
        copy(parents.begin(), parents.end(), pop.begin());
        // Copy from crossover_point...end to res.
        copy(mutants.begin(), mutants.end(), pop.begin() + parents.size());
    }

    vector<int> fit_last = fitness(weight, value, pop, capacity);
    cout << "Fitness da última geração: \n";
    print_vec(fit_last);
    int max_idx_last = idx_max(fit_last);
    return pop[max_idx_last];
}

// Generate n random integers between lower_range and upper_range.
void n_random_ints(vector<int> &int_vector, int n, int lower_range, int upper_range)
{
    // Create the random number with a diferent seed every time.
    myclock::duration d = myclock::now() - beginning;
    unsigned seed = d.count();
    default_random_engine generator(seed);
    uniform_int_distribution<int> distribution(lower_range, upper_range);

    int k;
    for (int i = 0; i < n; i++)
    {
        k = distribution(generator);
        int_vector.push_back(k);
    }
}

int main(int argc, char **argv)
{
    // Usage ./genetic_algo_knapsack [RANDOM: 0 or 1]
    int random = 1;
    int n_generations = 1000;
    int sol_per_pop = 16;
    for (int i = 1; i < argc; i++)
    {
        if (string(argv[i]) == "-r") {
            if (i + 1 < argc)
                random = stoi(argv[i + 1]);
            else
                cout << "Missing random value. Using default value: 1" << endl;
        } else if (string(argv[i]) == "-g")
        {
            if (i + 1 < argc)
                n_generations = stoi(argv[i + 1]);
            else
                cout << "Missing number of generations. Using default value: 1000" << endl;
        }
        if (string(argv[i]) == "-s")
        {
            if (i + 1 < argc)
                sol_per_pop = stoi(argv[i + 1]);
            else
                cout << "Missing number of solutions per population. Using default value: 16" << endl;
        }
        if (string(argv[i]) == "-h")
        {
            cout << "Usage: ./genetic_algo_knapsack [-r (0 | 1)] [-g (1..+inf)] [-s (1..+inf)]" << endl;
            return 0;
        }
    }

    int n_items, capacity;
    cout << "Enter the number of itens in knapsack: ";
    cin >> n_items;
    cout << "Enter the capacity of knapsack: ";
    cin >> capacity;
    cout << endl;

    // Allocating the knapsack.
    vector<int> weights, values, item_ids;
    item_ids.reserve(n_items);
    weights.reserve(n_items);
    values.reserve(n_items);

    // Initialize the vectors - Knapsack.
    if (random)
    {
        for (int i = 1; i < n_items; i++)
            item_ids.push_back(i);
        n_random_ints(weights, n_items, 1, 15);
        n_random_ints(values, n_items, 1, 100);
    }
    else
    {
        int w, v;
        for (int i = 0; i < n_items; i++)
        {
            cin >> w >> v;
            item_ids.push_back(i + 1);
            weights.push_back(w);
            values.push_back(v);
        }
    }

    // ----- Print Items -----
    cout << "Item\tPeso\tValor" << endl;
    for (int i = 0; i < n_items; i++)
        cout << item_ids[i] << '\t' << weights[i] << '\t' << values[i] << endl;

    // Generating initial population
    population initial_pop;
    initial_pop.reserve(sol_per_pop);

    // Initialize and print the population.
    for (int i = 0; i < sol_per_pop; i++)
    {
        vector<int> items_gene;
        items_gene.reserve(n_items);
        n_random_ints(items_gene, n_items, 0, 1);
        initial_pop.push_back(items_gene);
    }
    cout << "Inicial population: \n";
    print_population(initial_pop);

    // Run the genetic algo.
    vector<int> parameters = optimize(weights, values, initial_pop, sol_per_pop, n_generations, capacity);

    int w_sum = 0;
    int v_sum = 0;
    cout << "Itens included in the knapsack: " << endl;
    for (int i = 0; i < (int) parameters.size(); i++)
    {
        if (parameters[i] != 0)
        {
            cout << item_ids[i] << ' ';
            w_sum += weights[i];
            v_sum += values[i];
        }
    }
    cout << endl;
    cout << "Total knapsack weight: " << w_sum << endl;
    cout << "Total knapsack value: " << v_sum << endl;
}