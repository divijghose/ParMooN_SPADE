# BLAS and LAPACK routines used in the code 

## Monte Carlo Realization generation
### SVD of Covariance matrix 
#### Code Snippet

```MKL_INT m1 = N_U, n = N_U, lda = N_U, ldu = N_U, ldvt = N_U, info;
double superb[std::min(N_U, N_U) - 1];

double *S = new double[N_U];
double *U = new double[N_U * N_U];
double *Vt = new double[N_U * N_U];

info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', m1, n, C, lda,S, U, ldu, Vt, ldvt, superb);

cout << endl
        << endl;

if (info > 0)
{
    printf("The algorithm computing SVD failed to converge.\n");
    exit(1);
}
cout << " SVD COMPUTED" << endl;
```

#### Purpose
 DGESVD computes the singular value decomposition (SVD) of a real
 ```N_U x N_U``` matrix ```C```, also computing the left and/or right singular
 vectors. The SVD is written

      C = U * S * Vt



 __Note that the routine returns V**T, not V.__

#### Input
```MATRIX_LAYOUT```
indicates whether the input and output matrices are stored in row major order or column major order, where:\
matrix_layout = ```LAPACK_ROW_MAJOR```, the matrices are stored in row major order.

```JOBU``` is CHARACTER*1\
Specifies options for computing all or part of the matrix\
U := ```'A'```:  all M columns of U are returned in array

```JOBVT``` is CHARACTER*1 \
Specifies options for computing all or part of the matrix\
VT := ```'A'```:  all N rows of V**T are returned in the array VT;

```M``` is INTEGER \
The number of rows of the input matrix ```C``` (Covariance Matrix). \
```m1 = N_U```.

```N``` is INTEGER \
The number of columns of the input matrix  ```C``` (Covariance Matrix). \
```n = N_U```.

A is DOUBLE PRECISION array, dimension (LDA,N), in our case ```C``` (Covariance Matrix)\
On entry, the M-by-N matrix ```C```. \
On exit, 
if ```JOBU != 'O'``` and ```JOBVT != 'O'```, the contents of ```C```
are destroyed.

```LDA``` is INTEGER \
The leading dimension of the array ```C```. \
LDA >= max(1,M)\
```lda = N_U```

```S``` is DOUBLE PRECISION array, dimension (min(M,N))= ```N_U```.  
The singular values of A, sorted so that S(i) >= S(i+1).

```U ``` is DOUBLE PRECISION array, dimension ```LDU x N_U``` \
If ```JOBU = 'A'```, U contains the M-by-M orthogonal matrix U;

```LDU``` is INTEGER \
The leading dimension of the array U.\
```ldu = N_U```

	

```Vt``` is DOUBLE PRECISION array, ```LDVT x N_U```
If JOBVT = 'A', ```Vt``` contains the ```N_U x N_U``` orthogonal matrix V**T

```LDVT``` is INTEGER \
The leading dimension of the array VT. \
```ldvt = N_U```

#### Output
```S[i], 0<=i<N_U``` \
```U, N_U x N_U``` \
```Vt, N_U x N_U```