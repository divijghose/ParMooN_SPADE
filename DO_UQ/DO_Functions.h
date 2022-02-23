
// Convert array from row major to column major format for cache optimization
//`width` refers to number of columns of matrix and `height` refers to number of rows of matrix.
double* RowtoColMaj(double* a,int height,int width){
    double* b = new double[width*height]();
    for(int i=0;i<height;i++){
        for(int j=0;j<width;j++){
                b[height*j+i]=a[width*i+j];


        }
    }
    return b;
}

// Convert array from column major to row major format for cache optimization
//`width` refers to number of columns of matrix and `height` refers to number of rows of matrix.
double* ColtoRowMaj(double* a,int height,int width){
    double* b = new double[width*height]();
    for(int i=0;i<height;i++){
        for(int j=0;j<width;j++){
                b[width*i+j]=a[height*j+i];

        }
    }
    return b;
}

//Given a vector and its dimensions, calculate its Covariance. Default input and output are in row storage format.
//`width` refers to number of columns of matrix and `height` refers to number of rows of matrix.
double* CalcCovarianceMatx(double* Vector, int height, int width){
    double* Cov = new double[width* width]();
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,width,width, height, (1.0/(height-1)), Vector,height,Vector,width,0.0,Cov,width);

    return Cov;
}

//Given a vector and its dimensions, calculate its Coskewness . Default input and output are in Row storage format.
//`width` refers to number of columns of matrix and `height` refers to number of rows of matrix.
double* CalcCoskewnessMatx(double* Vector, int height, int width){

    double* M = new double[width*width*width]();
    for(int i = 0; i < width; i++){
        for(int j = 0; j < width; j++){
            for(int k = 0; k < width; k++){
                for(int p = 0; p < height; p++){

                    M[width*width*k + width*i + j] += ((Vector[width*p+i]*Vector[width*p+i]*Vector[width*p+i])/(height-1));

                }
            }
        }
    }

    return M;
}

lapack_int matInv(double *A, unsigned n)
{
    int ipiv[n+1];
    lapack_int ret;

    ret =  LAPACKE_dgetrf(LAPACK_COL_MAJOR,
                          n,
                          n,
                          A,
                          n,
                          ipiv);

    if (ret !=0)
        return ret;


    ret = LAPACKE_dgetri(LAPACK_COL_MAJOR,
                       n,
                       A,
                       n,
                       ipiv);
    return ret;
}

double* InvertCov(double* Cov, int N){
    double* CovInv = new double[N*N]();

    memcpy(CovInv,Cov, N*N*SizeOfDouble);
    matInv(CovInv,N);
    return CovInv;
    
} //change this to read Cov from database

