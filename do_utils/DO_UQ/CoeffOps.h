/**
 * @brief Routine to calculate the covariance matrix of coefficients and store it in TDatabase::ParamDB->COVARIANCE_MATRIX_DO
 * \f{equation}{
  C = \frac{1}{N_{R}-1}\cdot\phi^{T}\phi
\f}
 *
 * @param Vector Coefficient Matrix
 * @param height Number of Realizations, N_R (number of rows of Vector)
 * @param width Dimension of stochastic subspace, N_S (number of columns of Vector)
 *
 * @remark TDatabase::ParamDB->COVARIANCE_MATRIX_DO is defined in Database.C
 * @remark TDatabase::ParamDB->COVARIANCE_MATRIX_DO is stored in Row Major form
 */
void CalcCovarianceMatx(double *Vector)
{

    int height = TDatabase::ParamDB->REALISATIONS;
    int width = TDatabase::ParamDB->N_Subspace_Dim;
    // TDatabase::ParamDB->COVARIANCE_MATRIX_DO
    double *k = new double[width * width]();
    double *phi = new double[width * height](); // Col to Row Major
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            phi[width * i + j] = Vector[i + height * j];
        }
    }

    const double divVal = (1.0 / (height - 1));
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, width, width, height, divVal, phi, width, phi, width, 0.0, TDatabase::ParamDB->COVARIANCE_MATRIX_DO, width);
}


/**
 * @brief
 *
 * @param Vector Coefficient Matrix
 * @param height Number of Realizations, N_R (number of rows of Vector)
 * @param width Dimension of stochastic subspace, N_S (number of columns of Vector)
 *
 * @remark TDatabase::ParamDB->COVARIANCE_MATRIX_DO is defined in Database.C
 * @remark TDatabase::ParamDB->COVARIANCE_MATRIX_DO is stored in Row Major form
 */
void CalcCoskewnessMatx(double *Vector)
{

    int height = TDatabase::ParamDB->REALISATIONS;
    int width = TDatabase::ParamDB->N_Subspace_Dim;
    double *phi = new double[width * height](); // Col to Row Major
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            phi[width * i + j] = Vector[height * j + i];
        }
    }

    // TDatabase::ParamDB->COSKEWNESS_MATRIX_DO = new double[width * width * width]();
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < width; j++)
        {
            for (int k = 0; k < width; k++)
            {
                for (int p = 0; p < height; p++)
                {

                    TDatabase::ParamDB->COSKEWNESS_MATRIX_DO[width * width * k + width * i + j] += ((phi[width * p + i] * phi[width * p + j] * phi[width * p + k]) / (height - 1));
                }
            }
        }
    }
}


lapack_int matInv(double *A, unsigned n)
{
    int ipiv[n + 1];
    lapack_int ret;

    ret = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,
                         n,
                         n,
                         A,
                         n,
                         ipiv);

    if (ret != 0)
        // cout<<"Failure with dgetrf" << endl;
        return ret;

    ret = LAPACKE_dgetri(LAPACK_ROW_MAJOR,
                         n,
                         A,
                         n,
                         ipiv);

    // cout << "Failure with dgetri" << endl;
    return ret;
}

void InvertCov()
{
    int N = TDatabase::ParamDB->N_Subspace_Dim;
    // TDatabase::ParamDB->COVARIANCE_INVERSE_DO = new double[N * N]();

    memcpy(TDatabase::ParamDB->COVARIANCE_INVERSE_DO, TDatabase::ParamDB->COVARIANCE_MATRIX_DO, N * N * SizeOfDouble);
    matInv(TDatabase::ParamDB->COVARIANCE_INVERSE_DO, N);
}