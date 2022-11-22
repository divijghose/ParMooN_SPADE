#include <Database.h>
#include <FESpace2D.h>
#include <mkl.h>


int calculateStochSubspaceDim(TFESpace2D *Scalar_FeSpace, double *RealizationVector)
{
    int N_Realisations = TDatabase::ParamDB->REALIZATIONS;
    int N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
    int i;
    double *MeanVector = new double[N_DOF*1]();
    for (int i = 0; i < N_DOF; ++i)
    {
        for (int j = 0; j < N_Realisations; ++j)
        {
            MeanVector[i] += (RealizationVector[i * N_Realisations + j] / N_Realisations);
        }
    }
    double *PerturbationVector = new double[N_DOF * N_Realisations](); // \hat{C}^{i}_{dof} = C^{i}_{dof} - \overline{C}_{dof}
    for (int i = 0; i < N_DOF; ++i)
    {
        for (int j = 0; j < N_Realisations; ++j)
        {
            PerturbationVector[i * N_Realisations + j] = RealizationVector[i * N_Realisations + j] - MeanVector[i];
        }
    }
     //================================================================================================
    /////////////////////////////DO - Initialization SVD//////////////////////////////////////////////
    //================================================================================================
    // Declare SVD parameters
    int minDim = std::min(N_DOF, N_Realisations);
    MKL_INT mDO = N_DOF, nDO = N_Realisations, ldaDO = N_Realisations, lduDO = minDim, ldvtDO = N_Realisations, infoDO;
    double superbDO[minDim - 1];

    double *PerturbationVectorCopy = new double[N_DOF * N_Realisations]();
    memcpy(PerturbationVectorCopy, PerturbationVector, N_DOF * N_Realisations * SizeOfDouble);

    double *Sg = new double[minDim];
    double *L = new double[N_DOF * minDim];
    double *Rt = new double[minDim * N_Realisations];

    infoDO = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'S', 'N', mDO, nDO, PerturbationVectorCopy, ldaDO,
                            Sg, L, lduDO, Rt, ldvtDO, superbDO);

    if (infoDO > 0)
    {
        printf("The algorithm computing SVD for DO failed to converge.\n");
        exit(1);
    }
    cout << " DO SVD COMPUTED " << endl;

    //////////////////////////////////////////// DO - SVD End///////////////////////////////

    ///////DO - Subspace dimension calculation /////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////
    cout << "Starting DO Subspace Calculation" << endl;
    cout << "Subspace Manual " << TDatabase::ParamDB->Subspace_Manual << endl;
    cout << "Max subspace dim" << TDatabase::ParamDB->Max_Subspace_Dim << endl;
    int s1 = 0;
    int subDim = 0;
    if (TDatabase::ParamDB->Subspace_Manual == 1)
        subDim = TDatabase::ParamDB->Max_Subspace_Dim;
    else if (TDatabase::ParamDB->Subspace_Manual == 0)
    {
        double SVPercent = TDatabase::ParamDB->SVPERCENT;
        double valDO = 0.0;
        double sumSingularValDO = 0.0;
        for (int i = 0; i < minDim; i++)
        {
            sumSingularValDO += Sg[i];
        }
        while (valDO / sumSingularValDO < SVPercent)
        {
            valDO += Sg[s1];
            s1++;
        }

        subDim = s1 + 1;
        if (subDim > TDatabase::ParamDB->Max_Subspace_Dim)
            subDim = TDatabase::ParamDB->Max_Subspace_Dim;
    }
    else
    {
        cout << "Please enter correct value of Subspace_Manual (0 or 1)" << endl;
        exit(0);
    }

    cout << " SUBSPACE DIMENSION : " << subDim << endl;
    ////////Subspace dimension calculated//////////////////
    TDatabase::ParamDB->N_Subspace_Dim = subDim;

    return subDim;

}

void InitializeDO(TFESpace2D *Scalar_FeSpace, double *RealizationVector, double *MeanVector, double *ModeVector, double *CoeffVector)
{
    int N_Realisations = TDatabase::ParamDB->REALIZATIONS;
    int N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
    int subDim = TDatabase::ParamDB->N_Subspace_Dim;
    int i;
    for (int i = 0; i < N_DOF; ++i)
    {
        for (int j = 0; j < N_Realisations; ++j)
        {
            MeanVector[i] += (RealizationVector[i * N_Realisations + j] / N_Realisations);
        }
    }
    printToTxt("Init/Mean.txt", MeanVector, N_DOF, 1, 'C');

    double *PerturbationVector = new double[N_DOF * N_Realisations](); // \hat{C}^{i}_{dof} = C^{i}_{dof} - \overline{C}_{dof}
    for (int i = 0; i < N_DOF; ++i)
    {
        for (int j = 0; j < N_Realisations; ++j)
        {
            PerturbationVector[i * N_Realisations + j] = RealizationVector[i * N_Realisations + j] - MeanVector[i];
        }
    }
    printToTxt("Init/PerturbationVector.txt", PerturbationVector, N_DOF, N_Realisations, 'R');

     //================================================================================================
    /////////////////////////////DO - Initialization SVD//////////////////////////////////////////////
    //================================================================================================
    // Declare SVD parameters
    int minDim = std::min(N_DOF, N_Realisations);
    MKL_INT mDO = N_DOF, nDO = N_Realisations, ldaDO = N_Realisations, lduDO = minDim, ldvtDO = N_Realisations, infoDO;
    double superbDO[minDim - 1];

    double *PerturbationVectorCopy = new double[N_DOF * N_Realisations]();
    memcpy(PerturbationVectorCopy, PerturbationVector, N_DOF * N_Realisations * SizeOfDouble);

    double *Sg = new double[minDim];
    double *L = new double[N_DOF * minDim];
    double *Rt = new double[minDim * N_Realisations];

    infoDO = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'S', 'N', mDO, nDO, PerturbationVectorCopy, ldaDO,
                            Sg, L, lduDO, Rt, ldvtDO, superbDO);

    if (infoDO > 0)
    {
        printf("The algorithm computing SVD for DO failed to converge.\n");
        exit(1);
    }
    cout << " DO SVD COMPUTED " << endl;

    //////////////////////////////////////////// DO - SVD End///////////////////////////////

    /////Projection Matrix///////////
    ////////////////////////////////

    /// Initialize Coefficient Matrix - First subDim columns of Projection Matrix ////////////////////

    ////////////Initialize Mode Vector - First subDim columns of Left Singular Vector//////////////////


    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = 0; j < subDim; j++)
        {
            ModeVector[j * N_DOF + i] = L[i * minDim + j]; // ModeVector in Col Major form
        }
    }

    TFEVectFunct2D *FEFVector_ModeInit = new TFEVectFunct2D(Scalar_FeSpace, (char *)"C_Mode_Init", (char *)"Mode Component", ModeVector, N_DOF, subDim); 

    double *normzdModeVector = new double[N_DOF * subDim]();
    double *tempModeVector = new double[N_DOF * subDim]();
    normalizeStochasticModes(Scalar_FeSpace, FEFVector_ModeInit, subDim, tempModeVector);
    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = 0; j < subDim; j++)
        {
            normzdModeVector[i * subDim + j] = tempModeVector[j * N_DOF + i];
        }
    }
    double *ProjectionVector = new double[N_Realisations * subDim]();

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N_Realisations, subDim, N_DOF, 1.0, PerturbationVector, N_Realisations, normzdModeVector, subDim, 0.0, ProjectionVector, subDim);
    for (int i = 0; i < N_Realisations; i++)
    {
        for (int j = 0; j < subDim; j++)
        {
            CoeffVector[j * N_Realisations + i] = ProjectionVector[i * subDim + j]; // CoeffVector in Col Major form
        }
    }

    memcpy(ModeVector, tempModeVector, N_DOF * subDim * SizeOfDouble);
    printToTxt("Init/Coeff.txt", CoeffVector, N_Realisations, subDim, 'C');
    printToTxt("Init/Mode.txt", ModeVector, N_DOF, subDim, 'C');
    

    double *IPMatxModeInit = new double[subDim * subDim]();
    double *IPMatxMeanInit = new double[1 * 1]();

    calcIPMatx(IPMatxModeInit, ModeVector, N_DOF, subDim, 'C');
    printToTxt("Init/IPMatxMode_Init.txt", IPMatxModeInit, subDim, subDim, 'R');

    calcIPMatx(IPMatxMeanInit, MeanVector, N_DOF, 1, 'C');
    printToTxt("Init/IPMatxMean_Init.txt", IPMatxMeanInit, 1, 1, 'R');

    int m = 0;
    

    delete[] PerturbationVector;
    delete[] PerturbationVectorCopy;
    delete[] Sg;
    delete[] L;
    delete[] Rt;
    delete[] ProjectionVector;

    ////////////////////////////////////////////DO - Initialization Ends//////////////////////////////////////
    ///////================================================================================//////////////////
    return;
}