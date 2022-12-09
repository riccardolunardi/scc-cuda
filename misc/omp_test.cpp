#include <stdio.h>

#include <omp.h>

int main() {

    const int maxNumThreads = omp_get_max_threads();

    printf("Maximum number of threads for this machine: %i\n", maxNumThreads);

    printf("Not yet started a parallel Section: the number of threads is %i\n", omp_get_num_threads());

    printf("Setting the maximum number of threads...\n");
    omp_set_num_threads(maxNumThreads);

    printf("Once again, not yet started a parallel Section: the number of threads is still %i\n", omp_get_num_threads());

    printf("Starting a parallel Section...\n");

#pragma omp parallel for 
    for (int i = 0; i < maxNumThreads; i++) {
        int tid = omp_get_thread_num();
        printf("This is thread %i announcing that the number of launched threads is %i\n", tid, omp_get_num_threads());
    }

}