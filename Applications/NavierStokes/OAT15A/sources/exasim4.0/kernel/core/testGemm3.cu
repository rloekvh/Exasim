/* GEMM is a General Matrix Multiply - a subroutine in the Basic Linear Algebra Subprograms library*/

/* Includes, system */
#include <stdio.h>
#include <stdlib.h>
//#include <string.h>
//#include <cstring>
#include <sys/time.h>

/* Includes, cuda */
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "helper_cuda.h"
#include "gpuGEMM.cu"
#include "gpuArrayGemm.cu"

//using namespace std;

#define BLOCK_SIZE 16


void print1darray(double* a, int m)
{    
    for (int i=0; i<m; i++)
        printf("%g   ", a[i]);          
    printf("\n");
}

void print2darray(double* a, int m, int n)
{
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++)
            printf("%g   ", a[j*m+i]);              
        printf("\n");
    }
    printf("\n");
}

void print3darray(double* a, int m, int n, int p)
{
    for (int k=0; k<p; k++) {
        for (int i=0; i<m; i++) {
            for (int j=0; j<n; j++)
                printf("%g   ", a[k*n*m+j*m+i]);                  
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
}

void printArray2D(double* a, int m, int n)
{
    int N = m*n;
    double *b = (double*) malloc (sizeof (double)*N);
    cudaMemcpy(b, a, N*sizeof(double), cudaMemcpyDeviceToHost);    
    print2darray(b, m, n);
    free(b);
}

/* ======================================================= */
/* Simple host implementation of a simple version of sgemm */
/* ======================================================= */
static void simple_dgemm(int M, int N, int K, double alpha, const double *A, const double *B,
                         double beta, double *C) {
  int i, j, k;
  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j){
      double prod = 0;
      for (k = 0; k < K; ++k){
	prod += A[k * M + i] * B[j * K + k];
      }
      C[j * M + i] = alpha * prod + beta * C[j * M + i];
    }
  }
}

/* ======================= */
/* dgemm from BLAS library */
/* ======================= */
extern "C"{
extern void dgemm_(char *, char * , 
		  int *, int *, int *,
		  double *, double *, int *,
		  double *, int *,
		   double *, double *, int *); };

double geterror(double *C1, double *C2, int N)
{
    double *h_C1, *h_C2;
    h_C1 = (double *)malloc(N * sizeof(double) );  
    h_C2 = (double *)malloc(N * sizeof(double) );  
    
    //cublasGetVector(N, sizeof(double), C1, 1, h_C1, 1); 
    //cublasGetVector(N, sizeof(double), C2, 1, h_C2, 1);
    cudaMemcpy(h_C1, C1, N*sizeof(double), cudaMemcpyDeviceToHost);  
    cudaMemcpy(h_C2, C2, N*sizeof(double), cudaMemcpyDeviceToHost);  

    double errors = 0.0;
    for (int i=0; i<2; i++)
        if (fabs(h_C1[i]-h_C2[i])>errors)
            errors = fabs(h_C1[i]-h_C2[i]);

    //printArray2D(C1, 10, 10);
    //print2darray(h_C1, 10, 10);
    //printArray2D(C2, 10, 10);    
    //print2darray(h_C2, 10, 10);
   
    free(h_C1);
    free(h_C2);

    return errors;
}

static void gpuGemm(cublasHandle_t handle, double *times, double *errors, double *C, double *A, double *B, 
        double alpha, double beta, int I, int J, int K, double *Ctmp, double *Cblas, double *A1, double *A2, double *A3, 
        int *index, int I1, int K1, int I2, int K2, int I3, int K3, int elemtype)
{   
    struct timeval tv1, tv2;
    int M = I*K;

    cudaDeviceSynchronize();
    gettimeofday(&tv1, NULL);  
    cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, I, J, K, &alpha, A, I, B, K, &beta, Cblas, I);
    cudaDeviceSynchronize();
    gettimeofday(&tv2, NULL);
    times[0] = (double)(tv2.tv_usec-tv1.tv_usec)/1000 + (double)(tv2.tv_sec -tv1.tv_sec )*1000;   

    cudaDeviceSynchronize();
    gettimeofday(&tv1, NULL);  
    gpuGemmV0(C, A, B, alpha, beta, I, J, K);
    cudaDeviceSynchronize();
    gettimeofday(&tv2, NULL);
    times[1] = (double)(tv2.tv_usec-tv1.tv_usec)/1000 + (double)(tv2.tv_sec -tv1.tv_sec )*1000;
    errors[1] = geterror(C, Cblas, I*J);
  
    cudaDeviceSynchronize();
    gettimeofday(&tv1, NULL);  
    gpuGemmV3(C, A, B, alpha, beta, I, J, K);
    cudaDeviceSynchronize();
    gettimeofday(&tv2, NULL);
    times[2] = (double)(tv2.tv_usec-tv1.tv_usec)/1000 + (double)(tv2.tv_sec -tv1.tv_sec )*1000;
    errors[2] = geterror(C, Cblas, I*J);

    if (I<=128) {
        cudaDeviceSynchronize();
        gettimeofday(&tv1, NULL);  
        gpuGemmV7(C, A, B, alpha, beta, I, J, K);
        cudaDeviceSynchronize();
        gettimeofday(&tv2, NULL);
        times[3] = (double)(tv2.tv_usec-tv1.tv_usec)/1000 + (double)(tv2.tv_sec -tv1.tv_sec )*1000;
        errors[3] = geterror(C, Cblas, I*J);
    }

    if (M<=1280) {
        cudaDeviceSynchronize();
        gettimeofday(&tv1, NULL);  
        gpuGemmV2(C, A, B, alpha, beta, I, J, K);
        cudaDeviceSynchronize();
        gettimeofday(&tv2, NULL);
        times[4] = (double)(tv2.tv_usec-tv1.tv_usec)/1000 + (double)(tv2.tv_sec -tv1.tv_sec )*1000;
        errors[4] = geterror(C, Cblas, I*J);

        cudaDeviceSynchronize();
        gettimeofday(&tv1, NULL);  
        gpuGemmV5(C, A, B, alpha, beta, I, J, K);
        cudaDeviceSynchronize();
        gettimeofday(&tv2, NULL);
        times[5] = (double)(tv2.tv_usec-tv1.tv_usec)/1000 + (double)(tv2.tv_sec -tv1.tv_sec )*1000;
        errors[5] = geterror(C, Cblas, I*J);

        cudaDeviceSynchronize();
        gettimeofday(&tv1, NULL);  
        gpuGemmV9(C, A, B, alpha, beta, I, J, K);
        cudaDeviceSynchronize();
        gettimeofday(&tv2, NULL);
        times[6] = (double)(tv2.tv_usec-tv1.tv_usec)/1000 + (double)(tv2.tv_sec -tv1.tv_sec )*1000;
        errors[6] = geterror(C, Cblas, I*J);
    }

    if (M<=256) {
        cudaDeviceSynchronize();
        gettimeofday(&tv1, NULL);  
        gpuGemmV1(C, A, B, alpha, beta, I, J, K);
        cudaDeviceSynchronize();
        gettimeofday(&tv2, NULL);
        times[7] = (double)(tv2.tv_usec-tv1.tv_usec)/1000 + (double)(tv2.tv_sec -tv1.tv_sec )*1000;
        errors[7] = geterror(C, Cblas, I*J);

        cudaDeviceSynchronize();
        gettimeofday(&tv1, NULL);  
        gpuGemmV4(C, A, B, alpha, beta, I, J, K);
        cudaDeviceSynchronize();
        gettimeofday(&tv2, NULL);
        times[8] = (double)(tv2.tv_usec-tv1.tv_usec)/1000 + (double)(tv2.tv_sec -tv1.tv_sec )*1000;
        errors[8] = geterror(C, Cblas, I*J);

        cudaDeviceSynchronize();
        gettimeofday(&tv1, NULL);  
        gpuGemmV8(C, A, B, alpha, beta, I, J, K);
        cudaDeviceSynchronize();
        gettimeofday(&tv2, NULL);
        times[9] = (double)(tv2.tv_usec-tv1.tv_usec)/1000 + (double)(tv2.tv_sec -tv1.tv_sec )*1000;
        errors[9] = geterror(C, Cblas, I*J);
    }

    if (elemtype==1) {
        cudaDeviceSynchronize();
        gettimeofday(&tv1, NULL);  
        gpuTensorGemmV1(C, A1, A2, A3, B, Ctmp, index, I1, K1, I2, K2, I3, K3, J);
        cudaDeviceSynchronize();
        gettimeofday(&tv2, NULL);
        times[10] = (double)(tv2.tv_usec-tv1.tv_usec)/1000 + (double)(tv2.tv_sec -tv1.tv_sec )*1000;
        errors[10] = geterror(C, Cblas, I*J);

        cudaDeviceSynchronize();
        gettimeofday(&tv1, NULL);  
        gpuTensorGemmV4(C, A1, A2, A3, B, Ctmp, index, I1, K1, I2, K2, I3, K3, J);
        cudaDeviceSynchronize();
        gettimeofday(&tv2, NULL);
        times[11] = (double)(tv2.tv_usec-tv1.tv_usec)/1000 + (double)(tv2.tv_sec -tv1.tv_sec )*1000;
        errors[11] = geterror(C, Cblas, I*J);

        cudaDeviceSynchronize();
        gettimeofday(&tv1, NULL);  
        gpuTensorGemmV6(C, A1, A2, A3, B, Ctmp, index, I1, K1, I2, K2, I3, K3, J);
        cudaDeviceSynchronize();
        gettimeofday(&tv2, NULL);
        times[12] = (double)(tv2.tv_usec-tv1.tv_usec)/1000 + (double)(tv2.tv_sec -tv1.tv_sec )*1000;
        errors[12] = geterror(C, Cblas, I*J);

        cudaDeviceSynchronize();
        gettimeofday(&tv1, NULL);  
        gpuTensorGemmV8(C, A1, A2, A3, B, Ctmp, index, I1, K1, I2, K2, I3, K3, J);
        cudaDeviceSynchronize();
        gettimeofday(&tv2, NULL);
        times[13] = (double)(tv2.tv_usec-tv1.tv_usec)/1000 + (double)(tv2.tv_sec -tv1.tv_sec )*1000;
        errors[13] = geterror(C, Cblas, I*J);
    }
}


/* ==== */
/* Main */
/* ==== */
int main(int argc, char **argv)
{
  cublasStatus_t status;
  double *h_A, *h_B, *h_C, *h_D, *h_E, *h_F, *h_A1, *h_A2, *h_A3;
  double *d_A, *d_B, *d_C, *d_D, *d_E, *d_F, *d_A1, *d_A2, *d_A3;
  double *h_Cts, *d_Cts, *d_Ctm, *h_Fts, *d_Fts, *d_Ftm;
  int i, K, M1, N1, M2, N2, M3, N3, S1, S2, S3, SA, SB, SC, SD, SE, SF, SI;
  int *index1, *index2;    
  cublasHandle_t handle;
  struct timeval tv1, tv2;
  status = cublasCreate(&handle);
  double alpha = 1.0f;
  double beta = 0.0f;
  double times2[20];
  double errors2[20];
  double times3[20];
  double errors3[20];
  for (i=0; i<20; i++) {
    times2[i] = 1e10;
    errors2[i] = 0.0;
    times3[i] = 1e10;
    errors3[i] = 0.0;
  }

  K = 1024*16;
  
  M1 = 3; M2 = 3; M3 = 3;
  N1 = 3; N2 = 3; N3 = 3;
  M1 = 4; M2 = 4; M3 = 4;
  N1 = 4; N2 = 4; N3 = 4;

  S1 = M1*N1;
  S2 = M2*N2;
  S3 = M3*N3;
  SA = S1*S2;
  SB = N1*N2*K;  
  SC = M1*M2*K;
  SD = S1*S2*S3;  
  SE = N1*N2*N3*K;  
  SF = M1*M2*M3*K;  
  SI = M2*N1*K;

  h_A1 = (double *)malloc(S1 * sizeof(double) );  
  h_A2 = (double *)malloc(S2 * sizeof(double) );  
  h_A3 = (double *)malloc(S3 * sizeof(double) );  
  h_A = (double *)malloc(SA * sizeof(double) );  
  h_B = (double *)malloc(SB * sizeof(double) );  
  h_C = (double *)malloc(SC * sizeof(double) );  
  h_D = (double *)malloc(SD * sizeof(double) );  
  h_E = (double *)malloc(SE * sizeof(double) );  
  h_F = (double *)malloc(SF * sizeof(double) );  
  h_Cts = (double *)malloc(SC * sizeof(double) );  
  h_Fts = (double *)malloc(SF * sizeof(double) );

  for (i = 0; i < S1; i++)
    h_A1[i] = rand() / (double)RAND_MAX;
  for (i = 0; i < S2; i++)
    h_A2[i] = rand() / (double)RAND_MAX;
  for (i = 0; i < S3; i++)
    h_A3[i] = rand() / (double)RAND_MAX;
  for (i = 0; i < SA; i++)
    h_A[i] = rand() / (double)RAND_MAX;
  for (i = 0; i < SB; i++)
    h_B[i] = rand() / (double)RAND_MAX;
  for (i = 0; i < SC; i++)
    h_C[i] = rand() / (double)RAND_MAX;
  for (i = 0; i < SD; i++)
    h_D[i] = rand() / (double)RAND_MAX;
  for (i = 0; i < SE; i++)
    h_E[i] = rand() / (double)RAND_MAX;
  for (i = 0; i < SF; i++)
    h_F[i] = rand() / (double)RAND_MAX;

  cudaMalloc((void **)&d_A1, S1 * sizeof(double));
  cudaMalloc((void **)&d_A2, S2 * sizeof(double));
  cudaMalloc((void **)&d_A3, S3 * sizeof(double));
  cudaMalloc((void **)&d_A, SA * sizeof(double));
  cudaMalloc((void **)&d_B, SB * sizeof(double));
  cudaMalloc((void **)&d_C, SC * sizeof(double));
  cudaMalloc((void **)&d_D, SD * sizeof(double));
  cudaMalloc((void **)&d_E, SE * sizeof(double));
  cudaMalloc((void **)&d_F, SF * sizeof(double));
  cudaMalloc((void **)&d_Cts, SC * sizeof(double));
  cudaMalloc((void **)&d_Ctm, SC * sizeof(double));
  cudaMalloc((void **)&d_Fts, SF * sizeof(double));
  cudaMalloc((void **)&d_Ftm, SF * sizeof(double));  
  cudaMalloc((void **)&index1, SI * sizeof(int));
  cudaMalloc((void **)&index2, SC * sizeof(int));

  status = cublasSetVector(S1, sizeof(h_A1[0]), h_A1, 1, d_A1, 1);
  status = cublasSetVector(S2, sizeof(h_A2[0]), h_A2, 1, d_A2, 1);
  status = cublasSetVector(S3, sizeof(h_A3[0]), h_A3, 1, d_A3, 1);
  status = cublasSetVector(SA, sizeof(h_A[0]), h_A, 1, d_A, 1);
  status = cublasSetVector(SB, sizeof(h_B[0]), h_B, 1, d_B, 1);
  status = cublasSetVector(SC, sizeof(h_C[0]), h_C, 1, d_C, 1);
  status = cublasSetVector(SD, sizeof(h_D[0]), h_D, 1, d_D, 1);
  status = cublasSetVector(SE, sizeof(h_E[0]), h_E, 1, d_E, 1);
  status = cublasSetVector(SF, sizeof(h_F[0]), h_F, 1, d_F, 1);

  gpuKron(d_A, d_A1, d_A2, M1, N1, M2, N2);  
  gpuKron(d_D, d_A, d_A3, M1*M2, N1*N2, M3, N3);    

  gpuGemm(handle, times2, errors2, d_C, d_A, d_B, alpha, beta, M1*M2, K, N1*N2, d_Ctm, d_Cts, d_A1, d_A2, d_A3, index1, M1, N1, M2, N2, 0, 0, 1);
  for (i=0; i<14; i++)
      printf("Implementation: %d; 2D Execution time (in millisec): %.6f;  maximum error: %g\n", i, times2[i], errors2[i]);
    
  printf("\n\n");

  gpuGemm(handle, times3, errors3, d_F, d_D, d_E, alpha, beta, M1*M2*M3, K, N1*N2*N3, d_Ftm, d_Fts, d_A1, d_A2, d_A3, index1, M1, N1, M2, N2, M3, N3, 1);
  for (i=0; i<14; i++)
      printf("Implementation: %d; 3D Execution time (in millisec): %.6f;  maximum error: %g\n", i, times3[i], errors3[i]);

  cudaFree(d_A1);
  cudaFree(d_A2);
  cudaFree(d_A3);
  cudaFree(d_A);
  cudaFree(d_B);
  cudaFree(d_C);
  cudaFree(d_D);
  cudaFree(d_E);
  cudaFree(d_F);
  cudaFree(d_Cts);
  cudaFree(d_Ctm);
  
  free(h_A1);
  free(h_A2);
  free(h_A3);
  free(h_A);
  free(h_B);
  free(h_C);
  free(h_D);
  free(h_E);
  free(h_F);
  free(h_Cts);

  /* Shutdown */
  status = cublasDestroy(handle);

  return(0);
}
