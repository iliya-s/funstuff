#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mkl.h"

int main()
{
    double *A, *B, *C;
    int m, n, k, i, j;
    double alpha, beta;
    time_t initial, final;

    initial = time(NULL);

    printf("\nThis computes the real matrix C=alpha*A*B+beta*C using \n"
            "Intel mkl function dgemm, where A, B, and C are matrices and \n"
            "alpha and beta are double precision scalars\n\n");

    m = 15000, k = 15000, n = 15000;

    printf("Initialize data for matrix multiplication C=A*B for matrix \n"
            "A (%ix%i) and matrix B (%ix%i)\n\n", m, k, k, n);

    alpha = 1.0, beta = 0.0;

    printf("Allocating memory for matrices\n");
    
    A = mkl_malloc(m * k * sizeof(double), 64);
    B = mkl_malloc(k * n * sizeof(double), 64);
    C = mkl_malloc(m * n * sizeof(double), 64);

    if (A == NULL || B == NULL || C == NULL)
    {
        printf("\n Error: Can't allocate memory for matrices. Aborting... \n\n");

        mkl_free(A);
        mkl_free(B);
        mkl_free(C);

        return 1;
    }

    printf("Initialize matrix data\n");

    for (i = 0; i < (k*n); i++)
    {
        A[i] = (double) (i + 1);
    }

    for (i = 0; i < (k*n); i++)
    {
        B[i] = (double) (-i-1);
    }

    for (i = 0; i < (m*n); i++)
    {
        C[i] = 0.0;
    }

    printf("\nCompute matrix product using Intel MKL dgemm function via CBLAS\n");

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, k, B, n, beta, C, n);

    printf("\nComputation complete. \n");

    printf("\nDeallocating memory \n");
    mkl_free(A);
    mkl_free(B);
    mkl_free(C);

    printf("\nProgram Complete\n");
    
    final = time(NULL);

    printf("\nTotal time: %ld\n", final - initial);
    
    return 0;

}
