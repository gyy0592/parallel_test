#include <omp.h>
#include <iostream>
#include <chrono>

#include "benchmark.h"
int main() {

//   Fast test on whether omp is working.
  // Define loop dimensions
  int N = 10;
  int M = 20;

  // Define loop variables
  int i, j;

  // Use collapse clause to combine two loops into one parallel region
  #pragma omp parallel for collapse(2)
  for (i = 0; i < N; i++) {
    for (j = 0; j < M; j++) {
      // Do some computation inside the loops
      std::cout << "Thread " << omp_get_thread_num() << " is executing iteration (" << i << "," << j << ")\n";
    }
  }


//   Starting testing nested parfor 
// Enable nested parallelism
omp_set_nested(1);

auto start = std::chrono::high_resolution_clock::now();
// Create two tasks with different numbers of threads
#pragma omp parallel num_threads(2)
{
    // Get the task id (0 or 1)
    int task_id = omp_get_thread_num();

    // Run the benchmark function on each task with different numbers of threads
    #pragma omp parallel num_threads(2)
    {
        // Get the thread id within the task
        int thread_id = omp_get_thread_num();

        // Print some information
        printf("Task %d, thread %d\n", task_id, thread_id);

        // Call the benchmark function with some arguments
        benchmark1(10e11);
    }
}
auto end = std::chrono::high_resolution_clock::now();
auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
std::cout << "Nested parallel for took " << duration << " seconds\n";
return 0;
}