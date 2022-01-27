# BLAS and LAPACK routines used in the code 

## Contents ## {#contents}
1. [Monte Carlo realization generation](#monte-carlo-realization-generation) <br>
a. [SVD of Covariance Matrix](#svd-of-covariance-matrix)<br>
b. [Multiplication of truncated left singular value matrix with standard deviation matrix](#multiplication-of-truncated-left-singular-value-matrix-with-standard-deviation-matrix)
2. [DO Initialization](#do-initialization) <br>
a. [SVD of perturbation matrix](#svd-of-perturbation-matrix) <br>
b. [Calculation of Projection Matrix](#calculation-of-projection-matrix)
3. [DO Solver](#do-solver) <br>


## Monte Carlo Realization generation ## {#monte-carlo-realization-generation}
[Back to Contents](#contents)

### SVD of Covariance matrix ### {#svd-of-covariance-matrix}
[Back to Contents](#contents)

#### Code Snippet

```
MKL_INT m1 = N_U, n = N_U, lda = N_U, ldu = N_U, ldvt = N_U, info;
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
 `N_U x N_U` matrix `C`, also computing the left and/or right singular
 vectors. The SVD is written

      C = U * S * Vt



 __Note that the routine returns V**T, not V.__

#### Input
`MATRIX_LAYOUT`
indicates whether the input and output matrices are stored in row major order or column major order, where: <br>
matrix_layout = `LAPACK_ROW_MAJOR`, the matrices are stored in row major order.

`JOBU` is CHARACTER*1 <br>
Specifies options for computing all or part of the matrix <br>
U := ``'A'``:  all M columns of U are returned in array

`JOBVT` is CHARACTER*1 <br>
Specifies options for computing all or part of the matrix <br>
VT := ``'A'``:  all N rows of V**T are returned in the array VT;

`M` is INTEGER <br>
The number of rows of the input matrix `C` (Covariance Matrix). <br>
`m1 = N_U`.

`N` is INTEGER <br>
The number of columns of the input matrix  `C` (Covariance Matrix). <br>
`n = N_U`.

A is DOUBLE PRECISION array, dimension (LDA,N), in our case `C` (Covariance Matrix)<br>
On entry, the M-by-N matrix `C`. <br>
On exit, 
if ``JOBU != 'O'`` and ``JOBVT != 'O'``, the contents of `C`
are destroyed.

`LDA` is INTEGER <br>
The leading dimension of the array `C`. <br>
LDA >= max(1,M) <br>
`lda = N_U`

`S` is DOUBLE PRECISION array, dimension (min(M,N))= `N_U`.  
The singular values of A, sorted so that S(i) >= S(i+1).

`U ` is DOUBLE PRECISION array, dimension `LDU x N_U` <br>
If ``JOBU = 'A'``, U contains the M-by-M orthogonal matrix U;

`LDU` is INTEGER <br>
The leading dimension of the array U.<br>
`ldu = N_U`

	

`Vt` is DOUBLE PRECISION array, `LDVT x N_U`
If JOBVT = 'A', `Vt` contains the `N_U x N_U` orthogonal matrix V**T

`LDVT` is INTEGER <br>
The leading dimension of the array VT. <br>
`ldvt = N_U`

#### Output
`S[i], 0<=i<N_U` <br>
`U, N_U x N_U` <br>
`Vt, N_U x N_U`

_____________________________________________________
### Multiplication of truncated left singular value matrix with standard deviation matrix ### {#multiplication-of-truncated-left-singular-value-matrix-with-standard-deviation-matrix}
[Back to Contents](#contents)

#### Purpose
The dgemm routine calculates the product of double precision matrices <br>
 `C = alpha*(A x B) + k*C`<br>
Here, we want <br>
` RealizationVector = Ut x Z `

#### Code snippet
`cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,m, n, k, alpha, A, k, B, n, beta, C, n);`

#### Input 
`CblasRowMajor` <br>
Indicates that the matrices are stored in row major order, with the elements of each row of the matrix stored contiguously.

`CblasNoTrans` <br>
Enumeration type
indicating that the matrices A
and B
should not be transposed or conjugate transposed before multiplication. 

m, n, k
Integers indicating the size of the matrices:

A
: m
rows by k
columns<br>
B
: k
rows by n
columns<br>
C
: m
rows by n
columns<br>
Here, `A = Ut`, `B=Z`, `C = RealizationVector`, hence <br>
`m = N_U` <br>
`n = N_Realisations` <br>
`k = modDim`

alpha
Real value used to scale the product of matrices A
and B.<br>
`alpha = 1.0`

A <br>
Array used to store matrix A (`Ut` in our case)

k <br>
Leading dimension of array A (`Ut`)
, or the number of elements between successive rows (for row major storage)
in memory. <br>
`k = modDim`

B <br>
Array used to store matrix B (`Z` in our case)

n <br>
Leading dimension of array B (`Z`)
, or the number of elements between successive rows (for row major storage)
in memory. <br>
`n = N_Realisations`

beta <br>
Real value used to scale matrix C (`RealizationVector`)
`beta = 0.0`

C <br>
Array used to store matrix C (`RealizationVector` in our case)

n
Leading dimension of array C
, or the number of elements between successive rows (for row major storage)
in memory.
`n = N_Realisations`
____________________________________________________________________

## DO Initialization ## {#do-initialization}
[Back to Contents](#contents)
### SVD of perturbation matrix ### {#svd-of-perturbation-matrix}
[Back to Contents](#contents)
1. [If N_Realisations > N_U](#if-number-of-realisations-is-greater-than-number-of-degrees-of-freedom)
2. [If N_U > N_Realisations](#if-number-of-degrees-of-freedom-is-greater-than-number-of-realisations)

#### If number of realisations is greater than number of degrees of freedom #### {#if-number-of-realisations-is-greater-than-number-of-degrees-of-freedom}
[Back to SVD of perturbation matrix](#svd-of-perturbation-matrix) <br>
`min(N_R,N_U) =  N_U`

```
int minDim = std::min(N_U,N_Realisations);
//Here, minDim = N_U
MKL_INT mDO = N_U, nDO = N_Realisations, ldaDO = N_Realisations, lduDO = minDim, ldvtDO = N_Realisations, infoDO;
double superbDO[minDim-1];

double* PerturbationVectorCopy = new double[N_U * N_Realisations]();
memcpy(PerturbationVectorCopy,PerturbationVector,N_U*N_Realisations*SizeOfDouble);

double* Sg = new double[minDim];
double* L = new double[N_U*minDim];
double* Rt = new double[minDim*N_Realisations];

infoDO = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'S', 'N', mDO, nDO, PerturbationVectorCopy, ldaDO,
					Sg, L, lduDO, Rt, ldvtDO, superbDO );


if( infoDO > 0 ) {
	printf( "The algorithm computing SVD for DO failed to converge.\n" );
	exit( 1 );
}
cout << " DO SVD COMPUTED " <<endl;
```

#### Purpose
 DGESVD computes the singular value decomposition (SVD) of a real
 `N_U x N_Realisations` matrix `PerturbationVector`, also computing the left and/or right singular
 vectors. The SVD is written

      PerturbationVector = L * Sg * Rt



 __Note that the routine returns R**T, not R.__

#### Input
`MATRIX_LAYOUT`
indicates whether the input and output matrices are stored in row major order or column major order, where:<br>
matrix_layout = `LAPACK_ROW_MAJOR`, the matrices are stored in row major order.

`JOBU` is CHARACTER*1<br>
Specifies options for computing all or part of the matrix <br>
L := ``'S'``:  the first `min(N_U,N_Realisation)` columns of U (the left singular matrix)<br>
__Here we are assuming__ `N_Realisations > N_U`, __hence U will have dimension__ `N_U x N_U`


`JOBVT` is CHARACTER*1 <br>
Specifies options for computing all or part of the matrix<br>
VT := ``'N'``: no rows of V**T (no right singular vectors) are computed.<br>

``M` is INTEGER <br>
The number of rows of the input matrix `PerturbationVectorCopy` (Copy of `PerturbationVector`). <br>
`mDO = N_U`.

`N` is INTEGER <br>
The number of columns of the input matrix  `PerturbationVectorCopy`(Copy of `PerturbationVector`). <br>
`nDO = N_Realisations`.

A is DOUBLE PRECISION array, dimension (LDA,N), in our case `PerturbationVectorCopy` (Copy of `PerturbationVector`)<br>
On entry, the M-by-N matrix `PerturbationVectorCopy`. <br>
On exit, 
if ``JOBU != 'O'`` and ``JOBVT != 'O'``, the contents of `PerturbationVectorCopy`
are destroyed. (Hence we are using `PerturbationVectorCopy` and not `PerturbationVector`)

`LDA` is INTEGER <br>
The leading dimension of the array `PerturbationVectorCopy`. <br>
LDA >= max(1,M)<br>
`lda = N_Realisations` <br>
**Assumption of `N_Realisations > N_U` also applies here**

`Sg` is DOUBLE PRECISION array, dimension (min(M,N))= `N_U`.  
The singular values of A, sorted so that S(i) >= S(i+1).

`U ` is DOUBLE PRECISION array, dimension `LDU x N_U` <br>      
If JOBU = 'S', U contains the first min(m,n) columns of U
(the left singular vectors, stored columnwise)
          
`LDU` is INTEGER <br>
The leading dimension of the array U.<br>
`ldu = N_U`



#### Output
`Sg[i], 0<=i<N_U` <br>
`L, N_U x N_U` <br>
________________________________________________________

#### If number of degrees of freedom is greater than number of realisations #### {#if-number-of-degrees-of-freedom-is-greater-than-number-of-realisations}
[Back to SVD of perturbation matrix](#svd-of-perturbation-matrix) <br>
`min(N_R,N_U) =  N_Realisations`

```
int minDim = std::min(N_U,N_Realisations);
//Here, minDim = N_Realisations
MKL_INT mDO = N_U, nDO = N_Realisations, ldaDO = N_Realisations, lduDO = minDim, ldvtDO = N_Realisations, infoDO;
double superbDO[minDim-1];

double* PerturbationVectorCopy = new double[N_U * N_Realisations]();
memcpy(PerturbationVectorCopy,PerturbationVector,N_U*N_Realisations*SizeOfDouble);

double* Sg = new double[minDim];
double* L = new double[N_U*minDim];
double* Rt = new double[minDim*N_Realisations];

infoDO = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'S', 'N', mDO, nDO, PerturbationVectorCopy, ldaDO,
					Sg, L, lduDO, Rt, ldvtDO, superbDO );


if( infoDO > 0 ) {
	printf( "The algorithm computing SVD for DO failed to converge.\n" );
	exit( 1 );
}
cout << " DO SVD COMPUTED " <<endl;
```

#### Purpose
 DGESVD computes the singular value decomposition (SVD) of a real
 `N_U x N_Realisations` matrix `PerturbationVector`, also computing the left and/or right singular
 vectors. The SVD is written

      PerturbationVector = L * Sg * Rt



 __Note that the routine returns R**T, not R.__

#### Input
`MATRIX_LAYOUT`
indicates whether the input and output matrices are stored in row major order or column major order, where:<br>
matrix_layout = `LAPACK_ROW_MAJOR`, the matrices are stored in row major order.

`JOBU` is CHARACTER*1<br>
Specifies options for computing all or part of the matrix <br>
L := ``'S'`` :  the first `min(N_U,N_Realisation)` columns of U (the left singular matrix)<br>
__Here we are assuming__ `N_U > N_Realisations`, __hence U will have dimension__ `N_U x N_Realisations`


`JOBVT` is CHARACTER*1 <br>
Specifies options for computing all or part of the matrix<br>
VT := ``'N'``: no rows of V**T (no right singular vectors) are computed.<br>

`M` is INTEGER <br>
The number of rows of the input matrix `PerturbationVectorCopy`
(Copy of `PerturbationVector`). <br>
`mDO = N_U`.

`N` is INTEGER <br>
The number of columns of the input matrix  `PerturbationVectorCopy` (Copy of `PerturbationVector`). <br>
`nDO = N_Realisations`.

A is DOUBLE PRECISION array, dimension (LDA,N), in our case `PerturbationVectorCopy` (Copy of `PerturbationVector`)<br>
On entry, the M-by-N matrix `PerturbationVectorCopy`. <br>
On exit, 
if ``JOBU != 'O'`` and ``JOBVT != 'O'``, the contents of `PerturbationVectorCopy`
are destroyed.(Hence we are using `PerturbationVectorCopy` and not `PerturbationVector`)

`LDA` is INTEGER <br>
The leading dimension of the array `PerturbationVectorCopy`. <br>
LDA >= max(1,M)<br>
`lda = N_Realisations` <br>


`Sg` is DOUBLE PRECISION array, dimension (min(M,N))= `N_Realisations`.  
The singular values of A, sorted so that S(i) >= S(i+1).

`U ` is DOUBLE PRECISION array, dimension `LDU x N_U` <br>      
If JOBU = 'S', U contains the first min(m,n) columns of U
(the left singular vectors, stored columnwise)
          
`LDU` is INTEGER <br>
The leading dimension of the array U.<br>
`ldu = N_Realisations`



#### Output
`Sg[i], 0<=i<N_Realisations` <br>
`L, N_U x N_Realisations` 

________________________________________________________

### Calculation of Projection Matrix ### {#calculation-of-projection-matrix}
[Back to Contents](#contents)
#### Purpose
The dgemm routine calculates the product of double precision matrices <br>
`C = alpha*(A x B) + beta*C`<br>
Here, we want <br>
` ProjectionVector = PerturbationVector^T x L `

#### Code snippet
```
cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,N_Realisations,minDim,N_U,1.0,PerturbationVector,N_Realisations,L,minDim,0.0,ProjectionVector,minDim);
```

#### Input 
`CblasRowMajor` <br>
Indicates that the matrices are stored in row major order, with the elements of each row of the matrix stored contiguously.

`CblasTrans` <br>
Enumeration type
indicating that the matrix `Perturbation Vector` should be transposed before multiplication. 

`CblasNoTrans` <br>
Enumeration type indicating that `L`
should not be transposed or conjugate transposed before multiplication. 

`m, n, k`
Integers indicating the size of the matrices:

A
: m
rows by k
columns<br>
B
: k
rows by n
columns<br>
C
: m
rows by n
columns<br>
Here, `A = PerturbationVector^T`, `B=L`, `C = ProjectionVector`, hence <br>
`m = N_Realisations` <br>
`n = minDim` <br>
`k = N_U` 

alpha
Real value used to scale the product of matrices A
and B.<br>
`alpha = 1.0`

A <br>
Array used to store matrix A (`PerturbationVector` transposed using CblasTrans in our case)

k <br>
Leading dimension of array A (`PerturbationVector`)
, or the number of elements between successive rows (for row major storage)
in memory. <br>
`k = N_Realisations`

B <br>
Array used to store matrix B (`L` in our case)

n <br>
Leading dimension of array B (`L`)
, or the number of elements between successive rows (for row major storage)
in memory. <br>
`n = minDim`

beta <br>
Real value used to scale matrix C (`ProjectionVector`)
`beta = 0.0`

C <br>
Array used to store matrix C (`ProjectionVector` in our case)

n
Leading dimension of array C
, or the number of elements between successive rows (for row major storage)
in memory.
`n = minDim`

_________________________________________________________________________________________

## DO Solver ## {#do-solver}
[Back to Contents](#contents)

### Covariance Matrix for Coefficients ### {#covariance-matrix-coefficients}
[Back to Contents](#contents)
#### Purpose
The dgemm routine calculates the product of double precision matrices <br>
`C = alpha*(A x B) + beta*C`<br>
Here, we want <br>
` Cov = (1/(N_Realisations-1))* CoeffVector^T x CoeffVector `

#### Code snippet
```
cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,N_Realisations,minDim,N_U,1.0,PerturbationVector,N_Realisations,L,minDim,0.0,ProjectionVector,minDim);
```

#### Input 
`CblasRowMajor` <br>
Indicates that the matrices are stored in row major order, with the elements of each row of the matrix stored contiguously.

`CblasTrans` <br>
Enumeration type
indicating that the matrix `CoeffVector` should be transposed before multiplication. 

`CblasNoTrans` <br>
Enumeration type indicating that `CoeffVector`
should not be transposed or conjugate transposed before multiplication. 

`m, n, k`
Integers indicating the size of the matrices:

A
: m
rows by k
columns<br>
B
: k
rows by n
columns<br>
C
: m
rows by n
columns<br>
Here, `A = CoeffVector^T`, `B=CoeffVector`, `C = Cov`, hence <br>
`m = subDim` (Subspace Dimesnion) <br>
`n = subDim` <br>
`k = N_Realisations` 

alpha
Real value used to scale the product of matrices A
and B.<br>
`alpha = 1.0/(N_realisations-1)`

A <br>
Array used to store matrix A (`CoeffVector` transposed using CblasTrans in our case)

k <br>
Leading dimension of array A (`CoeffVector`^T)
, or the number of elements between successive rows (for row major storage)
in memory. <br>
`k = N_Realisations`

B <br>
Array used to store matrix B (`CoeffVector` in our case)

n <br>
Leading dimension of array B (`CoeffVector`)
, or the number of elements between successive rows (for row major storage)
in memory. <br>
`n = subDim`

beta <br>
Real value used to scale matrix C (`Cov`)
`beta = 0.0`

C <br>
Array used to store matrix C (`Cov` in our case)

n
Leading dimension of array C (`Cov` in our case)
, or the number of elements between successive rows (for row major storage)
in memory.
`n = subDim`

_________________________________________________________________________________________
****