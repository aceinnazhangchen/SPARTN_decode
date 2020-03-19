#ifndef _GNSS_MATH_H_
#define _GNSS_MATH_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

typedef int8_t bool_t;

#define matrix_at(A, x, y, rownum)    ((A)[(x) + (rownum) * (y)])

#define MI(i, j, n) ((i) * (n) + (j))

#define SMD(i) ((i)*((i)+1)/2)

#define SMI(i, j)	((i) > (j) ? (SMD(i) + (j)) : (SMD(j) + (i)))

#undef	ABS
#define ABS(a)	   (((a) < 0) ? -(a) : (a))

// set square matrix A to identity matrix with the rows of n
void eye(double *A, uint32_t n);
void matcpy(double *A, const double *B, uint32_t n, uint32_t m);

// matrix multiplication: C = alpha .*(A * B) + beta .* C 
void matmul(const char *tr, uint32_t n, uint32_t k, uint32_t m, double alpha,
	   const double *A, const double *B, double beta, double *C);

// matrix add: C = alpha * A + beta * B 
void matadd(double alpha, const double *A, double beta, const double *B, uint32_t n, uint32_t m, double *C);
void matminus_fast(const double *A, const double *B, uint32_t n, uint32_t m, double *C);
double norm(const double* a, int n);
double* mat(int n, int m);
int *imat(int n, int m);
double dot(const double * a, const double * b, int n);
int inv4(const double *a, double *b);

double mean_dat(double* dat, double n);
double median_dat(double* dat, int n);

int normv3(const double* a, double* b);
void cross3(const double* a, const double* b, double* c);


#ifdef __cplusplus
}
#endif
#endif