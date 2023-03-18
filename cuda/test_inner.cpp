#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cblas.h>

double dot_self_written(double *x, double *y, int n) {
    double res = 0;
    for (int i = 0; i < n; i++) {
        res += x[i] * y[i];
    }
    return res;
}

int main() {
    int n = 4000;
    double *x = malloc(n * sizeof(double));
    double *y = malloc(n * sizeof(double));

    // initialize arrays with random values
    for (int i = 0; i < n; i++) {
        x[i] = (double)rand() / RAND_MAX;
        y[i] = (double)rand() / RAND_MAX;
    }

    // time self-written inner product
    clock_t start = clock();
    double res1 = dot_self_written(x, y, n);
    clock_t end = clock();
    double time1 = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time for self-written inner product: %f seconds\n", time1);

    // time cblas_ddot
    start = clock();
    double res2 = cblas_ddot(n, x, 1, y, 1);
    end = clock();
    double time2 = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time for cblas_ddot: %f seconds\n", time2);

    printf("Results: %f (self-written) vs %f (cblas_ddot)\n", res1, res2);

    free(x);
    free(y);
    return 0;
}
