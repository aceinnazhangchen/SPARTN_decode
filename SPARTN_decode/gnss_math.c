#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "gnss_math.h"

#ifndef	FALSE
	#define	FALSE	(0)
#endif

#ifndef	TRUE
	#define	TRUE	(!FALSE)
#endif

#undef	MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))

#undef	MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))


#ifndef SIGN
#define SIGN(a)    ((a > 0) ? 1 : ((a < 0) ? -1 : 0))
#endif

#undef SQR
#define SQR(a)      ((a)*(a))

void
eye(double *A, uint32_t n)
{
    uint32_t i;
    
    for (i = 0; i < n; i++) 
		A[i + i * n] = 1.0;
}

static void *nav_memcpy(void *dest, const void *src, uint32_t count)
{
	int8_t *dst8 = (int8_t *)dest;
	int8_t *src8 = (int8_t *)src;

	if (count == 0 || dst8 == src8) {
		return dest;
	}

	while (count--) {
		*dst8++ = *src8++;
	}

	return dest;
}

extern void 
matcpy(double *A, const double *B, uint32_t n, uint32_t m)
{
    nav_memcpy(A, B, sizeof(double) * n * m);
}

/* multiply matrix -----------------------------------------------------------*/
void matmul(const char *tr, uint32_t n, uint32_t k, uint32_t m, double alpha,
	const double *A, const double *B, double beta, double *C)
{
	double d;
	uint32_t i, j, x;
	char f = tr[0] == 'N' ? (tr[1] == 'N' ? 1 : 2) : (tr[1] == 'N' ? 3 : 4);

	for (i = 0; i < n; i++) {
		for (j = 0; j < k; j++) {
			d = 0.0;
			switch (f) {
			case 1:
				for (x = 0; x < m; x++)
					d += A[i + x*n] * B[x + j*m];
				break;
			case 2:
				for (x = 0; x < m; x++)
					d += A[i + x*n] * B[j + x*k];
				break;
			case 3:
				for (x = 0; x < m; x++)
					d += A[x + i*m] * B[x + j*m];
				break;
			case 4:
				for (x = 0; x < m; x++)
					d += A[x + i*m] * B[j + x*k];
				break;
			}

			if (beta == 0.0)
				C[i + j*n] = alpha * d;
			else
				C[i + j*n] = alpha * d + beta * C[i + j*n];
		}
	}
}


void matadd(double alpha, const double *A, double beta, const double *B, uint32_t n, uint32_t m, double *C)
{
	uint32_t i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			C[i + j*n] = alpha * A[i + j*n] + beta * B[i + j*n];
		}
	}

}

void matminus_fast(const double *A, const double *B, uint32_t n, uint32_t m, double *C)
{
	uint32_t i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			C[i + j * n] = A[i + j * n] - B[i + j * n];
		}
	}
}

/* new integer matrix ----------------------------------------------------------
* allocate memory of integer matrix
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
int *imat(int n, int m)
{
    int *p;

    if (n <= 0 || m <= 0) return NULL;
    if (!(p = (int *)malloc(sizeof(int)*n*m))) {
        //  fatalerr("integer matrix memory allocation error: n=%d,m=%d\n",n,m);
    }
    return p;
}



extern double median_dat(double* dat, int n)
{
	int		 i, j;
	double	 temp;

	if (n == 0 || dat == NULL) return 0.0;
	/* sort */
	for (i = 0; i < n; ++i) {
		for (j = i + 1; j < n; ++j) {
			if (dat[i] > dat[j]) {
				temp = dat[i];
				dat[i] = dat[j];
				dat[j] = temp;
			}
		}
	}
	if (n % 2 == 0) {
		return (dat[n / 2 - 1] + dat[n / 2]) / 2.0;
	}
	else {
		return dat[n / 2];
	}
}

extern double mean_dat(double* dat, double n)
{
	int		 i;
	double	 sum = 0.0;

	if (n == 0 || dat == NULL) return 0.0;
	/* sort */
	for (i = 0; i < n; ++i)
	{
		sum += dat[i];
	}

	return sum / n;
}

#define MATRIX_EPSILON (1e-60)
extern int inv4(const double* a, double* b) {
	double det =
		(((a[4 * 1 + 0] * -((a[4 * 2 + 1] * -(a[4 * 0 + 2] * a[4 * 3 + 3] -
			a[4 * 0 + 3] * a[4 * 3 + 2]) +
			a[4 * 2 + 2] * (a[4 * 0 + 1] * a[4 * 3 + 3] -
				a[4 * 0 + 3] * a[4 * 3 + 1])) +
			a[4 * 2 + 3] * -(a[4 * 0 + 1] * a[4 * 3 + 2] -
				a[4 * 0 + 2] * a[4 * 3 + 1])) +
			a[4 * 1 + 1] * ((a[4 * 2 + 0] * -(a[4 * 0 + 2] * a[4 * 3 + 3] -
				a[4 * 0 + 3] * a[4 * 3 + 2]) +
				a[4 * 2 + 2] * (a[4 * 0 + 0] * a[4 * 3 + 3] -
					a[4 * 0 + 3] * a[4 * 3 + 0])) +
				a[4 * 2 + 3] * -(a[4 * 0 + 0] * a[4 * 3 + 2] -
					a[4 * 0 + 2] * a[4 * 3 + 0]))) +
			a[4 * 1 + 2] * -((a[4 * 2 + 0] * -(a[4 * 0 + 1] * a[4 * 3 + 3] -
				a[4 * 0 + 3] * a[4 * 3 + 1]) +
				a[4 * 2 + 1] * (a[4 * 0 + 0] * a[4 * 3 + 3] -
					a[4 * 0 + 3] * a[4 * 3 + 0])) +
				a[4 * 2 + 3] * -(a[4 * 0 + 0] * a[4 * 3 + 1] -
					a[4 * 0 + 1] * a[4 * 3 + 0]))) +
			a[4 * 1 + 3] * ((a[4 * 2 + 0] * -(a[4 * 0 + 1] * a[4 * 3 + 2] -
				a[4 * 0 + 2] * a[4 * 3 + 1]) +
				a[4 * 2 + 1] * (a[4 * 0 + 0] * a[4 * 3 + 2] -
					a[4 * 0 + 2] * a[4 * 3 + 0])) +
				a[4 * 2 + 2] * -(a[4 * 0 + 0] * a[4 * 3 + 1] -
					a[4 * 0 + 1] * a[4 * 3 + 0])));

	if (fabs(det) < MATRIX_EPSILON) return -1;

	b[4 * 0 + 0] =
		((a[4 * 2 + 1] *
			-(a[4 * 1 + 2] * a[4 * 3 + 3] - a[4 * 1 + 3] * a[4 * 3 + 2]) +
			a[4 * 2 + 2] *
			(a[4 * 1 + 1] * a[4 * 3 + 3] - a[4 * 1 + 3] * a[4 * 3 + 1])) +
			a[4 * 2 + 3] *
			-(a[4 * 1 + 1] * a[4 * 3 + 2] - a[4 * 1 + 2] * a[4 * 3 + 1])) /
		det;
	b[4 * 1 + 0] =
		-((a[4 * 2 + 0] *
			-(a[4 * 1 + 2] * a[4 * 3 + 3] - a[4 * 1 + 3] * a[4 * 3 + 2]) +
			a[4 * 2 + 2] *
			(a[4 * 1 + 0] * a[4 * 3 + 3] - a[4 * 1 + 3] * a[4 * 3 + 0])) +
			a[4 * 2 + 3] *
			-(a[4 * 1 + 0] * a[4 * 3 + 2] - a[4 * 1 + 2] * a[4 * 3 + 0])) /
		det;
	b[4 * 2 + 0] =
		((a[4 * 2 + 0] *
			-(a[4 * 1 + 1] * a[4 * 3 + 3] - a[4 * 1 + 3] * a[4 * 3 + 1]) +
			a[4 * 2 + 1] *
			(a[4 * 1 + 0] * a[4 * 3 + 3] - a[4 * 1 + 3] * a[4 * 3 + 0])) +
			a[4 * 2 + 3] *
			-(a[4 * 1 + 0] * a[4 * 3 + 1] - a[4 * 1 + 1] * a[4 * 3 + 0])) /
		det;
	b[4 * 3 + 0] =
		-((a[4 * 2 + 0] *
			-(a[4 * 1 + 1] * a[4 * 3 + 2] - a[4 * 1 + 2] * a[4 * 3 + 1]) +
			a[4 * 2 + 1] *
			(a[4 * 1 + 0] * a[4 * 3 + 2] - a[4 * 1 + 2] * a[4 * 3 + 0])) +
			a[4 * 2 + 2] *
			-(a[4 * 1 + 0] * a[4 * 3 + 1] - a[4 * 1 + 1] * a[4 * 3 + 0])) /
		det;

	b[4 * 0 + 1] =
		-((a[4 * 2 + 1] *
			-(a[4 * 0 + 2] * a[4 * 3 + 3] - a[4 * 0 + 3] * a[4 * 3 + 2]) +
			a[4 * 2 + 2] *
			(a[4 * 0 + 1] * a[4 * 3 + 3] - a[4 * 0 + 3] * a[4 * 3 + 1])) +
			a[4 * 2 + 3] *
			-(a[4 * 0 + 1] * a[4 * 3 + 2] - a[4 * 0 + 2] * a[4 * 3 + 1])) /
		det;
	b[4 * 1 + 1] =
		((a[4 * 2 + 0] *
			-(a[4 * 0 + 2] * a[4 * 3 + 3] - a[4 * 0 + 3] * a[4 * 3 + 2]) +
			a[4 * 2 + 2] *
			(a[4 * 0 + 0] * a[4 * 3 + 3] - a[4 * 0 + 3] * a[4 * 3 + 0])) +
			a[4 * 2 + 3] *
			-(a[4 * 0 + 0] * a[4 * 3 + 2] - a[4 * 0 + 2] * a[4 * 3 + 0])) /
		det;
	b[4 * 2 + 1] =
		-((a[4 * 2 + 0] *
			-(a[4 * 0 + 1] * a[4 * 3 + 3] - a[4 * 0 + 3] * a[4 * 3 + 1]) +
			a[4 * 2 + 1] *
			(a[4 * 0 + 0] * a[4 * 3 + 3] - a[4 * 0 + 3] * a[4 * 3 + 0])) +
			a[4 * 2 + 3] *
			-(a[4 * 0 + 0] * a[4 * 3 + 1] - a[4 * 0 + 1] * a[4 * 3 + 0])) /
		det;
	b[4 * 3 + 1] =
		((a[4 * 2 + 0] *
			-(a[4 * 0 + 1] * a[4 * 3 + 2] - a[4 * 0 + 2] * a[4 * 3 + 1]) +
			a[4 * 2 + 1] *
			(a[4 * 0 + 0] * a[4 * 3 + 2] - a[4 * 0 + 2] * a[4 * 3 + 0])) +
			a[4 * 2 + 2] *
			-(a[4 * 0 + 0] * a[4 * 3 + 1] - a[4 * 0 + 1] * a[4 * 3 + 0])) /
		det;

	b[4 * 0 + 2] =
		((a[4 * 1 + 1] *
			-(a[4 * 0 + 2] * a[4 * 3 + 3] - a[4 * 0 + 3] * a[4 * 3 + 2]) +
			a[4 * 1 + 2] *
			(a[4 * 0 + 1] * a[4 * 3 + 3] - a[4 * 0 + 3] * a[4 * 3 + 1])) +
			a[4 * 1 + 3] *
			-(a[4 * 0 + 1] * a[4 * 3 + 2] - a[4 * 0 + 2] * a[4 * 3 + 1])) /
		det;
	b[4 * 1 + 2] =
		-((a[4 * 1 + 0] *
			-(a[4 * 0 + 2] * a[4 * 3 + 3] - a[4 * 0 + 3] * a[4 * 3 + 2]) +
			a[4 * 1 + 2] *
			(a[4 * 0 + 0] * a[4 * 3 + 3] - a[4 * 0 + 3] * a[4 * 3 + 0])) +
			a[4 * 1 + 3] *
			-(a[4 * 0 + 0] * a[4 * 3 + 2] - a[4 * 0 + 2] * a[4 * 3 + 0])) /
		det;
	b[4 * 2 + 2] =
		((a[4 * 1 + 0] *
			-(a[4 * 0 + 1] * a[4 * 3 + 3] - a[4 * 0 + 3] * a[4 * 3 + 1]) +
			a[4 * 1 + 1] *
			(a[4 * 0 + 0] * a[4 * 3 + 3] - a[4 * 0 + 3] * a[4 * 3 + 0])) +
			a[4 * 1 + 3] *
			-(a[4 * 0 + 0] * a[4 * 3 + 1] - a[4 * 0 + 1] * a[4 * 3 + 0])) /
		det;
	b[4 * 3 + 2] =
		-((a[4 * 1 + 0] *
			-(a[4 * 0 + 1] * a[4 * 3 + 2] - a[4 * 0 + 2] * a[4 * 3 + 1]) +
			a[4 * 1 + 1] *
			(a[4 * 0 + 0] * a[4 * 3 + 2] - a[4 * 0 + 2] * a[4 * 3 + 0])) +
			a[4 * 1 + 2] *
			-(a[4 * 0 + 0] * a[4 * 3 + 1] - a[4 * 0 + 1] * a[4 * 3 + 0])) /
		det;

	b[4 * 0 + 3] =
		-((a[4 * 1 + 1] *
			-(a[4 * 0 + 2] * a[4 * 2 + 3] - a[4 * 0 + 3] * a[4 * 2 + 2]) +
			a[4 * 1 + 2] *
			(a[4 * 0 + 1] * a[4 * 2 + 3] - a[4 * 0 + 3] * a[4 * 2 + 1])) +
			a[4 * 1 + 3] *
			-(a[4 * 0 + 1] * a[4 * 2 + 2] - a[4 * 0 + 2] * a[4 * 2 + 1])) /
		det;
	b[4 * 1 + 3] =
		((a[4 * 1 + 0] *
			-(a[4 * 0 + 2] * a[4 * 2 + 3] - a[4 * 0 + 3] * a[4 * 2 + 2]) +
			a[4 * 1 + 2] *
			(a[4 * 0 + 0] * a[4 * 2 + 3] - a[4 * 0 + 3] * a[4 * 2 + 0])) +
			a[4 * 1 + 3] *
			-(a[4 * 0 + 0] * a[4 * 2 + 2] - a[4 * 0 + 2] * a[4 * 2 + 0])) /
		det;
	b[4 * 2 + 3] =
		-((a[4 * 1 + 0] *
			-(a[4 * 0 + 1] * a[4 * 2 + 3] - a[4 * 0 + 3] * a[4 * 2 + 1]) +
			a[4 * 1 + 1] *
			(a[4 * 0 + 0] * a[4 * 2 + 3] - a[4 * 0 + 3] * a[4 * 2 + 0])) +
			a[4 * 1 + 3] *
			-(a[4 * 0 + 0] * a[4 * 2 + 1] - a[4 * 0 + 1] * a[4 * 2 + 0])) /
		det;
	b[4 * 3 + 3] =
		((a[4 * 1 + 0] *
			-(a[4 * 0 + 1] * a[4 * 2 + 2] - a[4 * 0 + 2] * a[4 * 2 + 1]) +
			a[4 * 1 + 1] *
			(a[4 * 0 + 0] * a[4 * 2 + 2] - a[4 * 0 + 2] * a[4 * 2 + 0])) +
			a[4 * 1 + 2] *
			-(a[4 * 0 + 0] * a[4 * 2 + 1] - a[4 * 0 + 1] * a[4 * 2 + 0])) /
		det;

	return 0;
}

/* new matrix ------------------------------------------------------------------
* allocate memory of matrix
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
double* mat(int n, int m)
{
	double* p;

	if (n <= 0 || m <= 0) return NULL;
	p = (double*)malloc(sizeof(double) * n * m);

	return p;
}


/* inner product ---------------------------------------------------------------
* inner product of vectors
* args   : double *a,*b     I   vector a,b (n x 1)
*          int    n         I   size of vector a,b
* return : a'*b
*-----------------------------------------------------------------------------*/
extern double dot(const double* a, const double* b, int n)
{
	double c = 0.0;

	while (--n >= 0) c += a[n] * b[n];
	return c;
}
/* euclid norm -----------------------------------------------------------------
* euclid norm of vector
* args   : double *a        I   vector a (n x 1)
*          int    n         I   size of vector a
* return : || a ||
*-----------------------------------------------------------------------------*/
extern double norm(const double* a, int n)
{
	return sqrt(dot(a, a, n));
}
/* outer product of 3d vectors -------------------------------------------------
* outer product of 3d vectors
* args   : double *a,*b     I   vector a,b (3 x 1)
*          double *c        O   outer product (a x b) (3 x 1)
* return : none
*-----------------------------------------------------------------------------*/
extern void cross3(const double* a, const double* b, double* c)
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}
/* normalize 3d vector ---------------------------------------------------------
* normalize 3d vector
* args   : double *a        I   vector a (3 x 1)
*          double *b        O   normlized vector (3 x 1) || b || = 1
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int normv3(const double* a, double* b)
{
	double r;
	if ((r = norm(a, 3)) <= 0.0) return 0;
	b[0] = a[0] / r;
	b[1] = a[1] / r;
	b[2] = a[2] / r;
	return 1;
}