/**
 * @file TLA2D_ParMooN_IVUQ_DO_Test.C
 * @brief Purpose:     Main program for scalar equations with new
 *                     kernels of ParMooN, for solution of
 *                     time-dependent linear advection equation
 *                     with uncertainty quantification
                       Features included in this main program -
                       1. Monte Carlo Realization Generation for scalar quantity of interest
                       2. Initialization of Mean, Modes and Coefficients for solution of Dynamically Orthogonal System of Equations
                       3. Solution of Mean and Mode PDEs and Coefficient ODEs.


 * @author Sashikumaar Ganesan
 * @author Divij Ghose
 * @author Thivin Anandh

 *
 */

// =======================================================================
//
// Purpose:     Main program for scalar equations with new kernels of ParMooN, for solution of time-dependent linear advection equation with
//              uncertainty quantification
//              Features included in this main program -
//              1. Monte Carlo Realization Generation for scalar quantity of interest
//              2. Initialization of Mean, Modes and Coefficients for solution of Dynamically Orthogonal System of Equations
//              3. Solution of Mean and Mode PDEs and Coefficient ODEs.
//
// Authors:      Sashikumaar Ganesan, Divij Ghose, Thivin Anandh
//
// History:     Implementation started on 10.03.2022

// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <SystemTCD2D.h>
#include <SystemTCD2D_ALE.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>
#include <TimeDiscRout.h>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <mkl.h>
#include <cmath>
#include <random>
#include <FEVectFunct2D.h>

// =======================================================================
// include current example
// =======================================================================
#include "../Examples/DO_UQ/linear_advection_do_test.h"

int main(int argc, char *argv[])
{
    int i, j, l, m, N_SubSteps, ORDER, N_Cells, N_DOF, img = 1, N_G;
    int N_Active;

    double *sol, *rhs, *oldrhs, t1, t2, errors[5], Linfty;
    double tau, end_time, *defect, olderror, olderror1, hmin, hmax;

    double *solMean, *rhsMean, *old_rhsMean;

    bool UpdateStiffnessMat, UpdateRhs, ConvectionFirstTime;

    char *VtkBaseName;
    char *VtkBaseNameMean;
    char *VtkBaseNameMode;

    const char vtkdir[] = "VTK";

    TDomain *Domain;
    TDatabase *Database = new TDatabase();
    TFEDatabase2D *FEDatabase = new TFEDatabase2D();
    TCollection *coll;
    TFESpace2D *Scalar_FeSpace, *fesp[1];
    TFEFunction2D *Scalar_FeFunction_Mean;
    TFEFunction2D *Scalar_FeFunction;
    TOutput2D *Output;
    TOutput2D *OutputMean;
    TOutput2D *OutputMode;
    TSystemTCD2D *SystemMatrix;
    TAuxParam2D *aux;
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};

    std::ostringstream os;
    os << " ";

    // ======================================================================
    // set the database values and generate mesh
    // ======================================================================
    // set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based
    Domain = new TDomain(argv[1]);

    if (TDatabase::ParamDB->PROBLEM_TYPE == 0)
        TDatabase::ParamDB->PROBLEM_TYPE = 2;
    OpenFiles();

    Database->WriteParamDB(argv[0]);
    Database->WriteTimeDB();
    ExampleFile();

    /* include the mesh from a mesh generator, for a standard mesh use the
     * build-in function. The GEOFILE describes the boundary of the domain. */
    if (TDatabase::ParamDB->MESH_TYPE == 0)
    {
        Domain->ReadGeo(TDatabase::ParamDB->GEOFILE);
        OutPut("PRM-GEO used for meshing !!!" << endl);
    } // ParMooN  build-in Geo mesh
    else if (TDatabase::ParamDB->MESH_TYPE == 1)
    {
        Domain->GmshGen(TDatabase::ParamDB->GEOFILE);
        OutPut("GMSH used for meshing !!!" << endl);
    }                                            // gmsh mesh
    else if (TDatabase::ParamDB->MESH_TYPE == 2) // ngle mesh
    {
        OutPut("Triangle.h used for meshing !!!" << endl);
        // TriaReMeshGen(Domain);
    }
    else
    {
        OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
        exit(0);
    }

#if defined(__HEMKER__) || defined(__BEAM__)
    TriaReMeshGen(Domain);
    TDatabase::ParamDB->UNIFORM_STEPS = 0;
#endif

    // refine grid up to the coarsest level
    for (i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
        Domain->RegRefineAll();

    // write grid into an Postscript file
    if (TDatabase::ParamDB->WRITE_PS)
        Domain->PS("Domain.ps", It_Finest, 0);

    // create output directory, if not already existing
    mkdir(vtkdir, 0777);
    const char modedir[] = "Modes";
    const char meandir[] = "Mean";
    const char coeffdir[] = "Coefficients";
    const char mcdir[] = "MonteCarlo";
    const char endir[] = "Energy_Data";

    mkdir(meandir, 0777);
    mkdir(modedir, 0777);
    mkdir(coeffdir, 0777);
    mkdir(mcdir, 0777);
    mkdir(endir, 0777);

    //=========================================================================
    // construct all finite element spaces
    //=========================================================================
    ORDER = TDatabase::ParamDB->ANSATZ_ORDER;

    coll = Domain->GetCollection(It_Finest, 0);
    N_Cells = coll->GetN_Cells();
    OutPut("N_Cells (space) : " << N_Cells << endl);

    // fespaces for scalar equation
    Scalar_FeSpace = new TFESpace2D(coll, (char *)"FE Space", (char *)"Solution Space",
                                    BoundCondition, 1, NULL);

    N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
    N_Active = Scalar_FeSpace->GetActiveBound();
    OutPut("dof all      : " << setw(10) << N_DOF << endl);
    OutPut("dof active   : " << setw(10) << N_Active << endl);

    // ///////////////////////////////////////////////////////////////////////////////////////////////
    // ////////// -------- REALISATION DATA GENERATION ----------------------------------------- //////
    // ///////////////////////////////////////////////////////////////////////////////////////////////

    int N_Realisations = TDatabase::ParamDB->REALIZATIONS;
    double LengthScale = TDatabase::ParamDB->LENGTHSCALE;
    double EigenPercent = TDatabase::ParamDB->EIGENPERCENT;

    double *org_x_coord = new double[N_DOF];
    double *org_y_coord = new double[N_DOF];
    double *x_coord = new double[N_DOF];
    double *y_coord = new double[N_DOF];
    int *mappingArray = new int[N_DOF];

    i = 0;
    int N = pow(2, TDatabase::ParamDB->UNIFORM_STEPS) + 1;
    for (int i = 0; i < N_DOF; i++)
    {
        int local_i = i / N;
        int local_j = i % N;

        x_coord[i] = double(1.0 / (N - 1)) * local_i;
        y_coord[i] = double(1.0 / (N - 1)) * local_j;
    }

    // cout << " End File Read" << endl;

    Scalar_FeSpace->GetDOFPosition(org_x_coord, org_y_coord);

    for (int i = 0; i < N_DOF; i++) // Generated Values
    {
        // get the generated Value
        double xx = x_coord[i];
        double yy = y_coord[i];
        bool foundFlag = false;

        for (int j = 0; j < N_DOF; j++) // Actual parmooN Co-ordinates
        {
            if (abs(xx - org_x_coord[j]) < 1e-10 && abs(yy - org_y_coord[j]) < 1e-10)
            {
                mappingArray[i] = j;
                foundFlag = true;
            }
        }

        if (!foundFlag)
            cerr << " DOF NOT FOUND FOR " << i << " position : " << setw(8) << org_x_coord[i] << setw(8) << org_y_coord[i] << endl;
    }
    // int N_DOF =  N * N;
    double *x = new double[N_DOF];
    double *y = new double[N_DOF];

    for (int i = 0; i < N_DOF; i++)
    {
        int local_i = i / N;
        int local_j = i % N;

        x[i] = double(1.0 / (N - 1)) * local_j;
        y[i] = double(1.0 / (N - 1)) * local_i;
    }

    double *C = new double[N_DOF * N_DOF];  // MATRIX
    double *C1 = new double[N_DOF * N_DOF]; // MATRIX  - Corelation Matrix

    double norm = 0;
    for (int i = 0; i < N_DOF; i++)
    {
        double actual_x = x[i];
        double actual_y = y[i];

        for (int j = 0; j < N_DOF; j++)
        {
            double local_x = x[j];
            double local_y = y[j];

            double r = sqrt(pow((actual_x - local_x), 2) + pow((actual_y - local_y), 2));

            // CO -Relation
            C[j * N_DOF + i] = exp((-1.0 * r) / (LengthScale))*(1+(r/LengthScale)+pow(r/LengthScale,2));
            C1[j * N_DOF + i] = exp((-1.0 * r) / (LengthScale))*(1+(r/LengthScale)+pow(r/LengthScale,2));

            if (TDatabase::ParamDB->stddev_switch == 0)
            {
                double sig_r1 = exp(-1.0 / (1.0 - pow((2 * actual_x - 1), 4))) * exp(-1.0 / (1 - pow((2 * actual_y - 1), 4)));
                double sig_r2 = exp(-1.0 / (1.0 - pow((2 * local_x - 1), 4))) * exp(-1.0 / (1 - pow((2 * local_y - 1), 4)));

                // Co Variance
                C[j * N_DOF + i] *= sig_r1 * sig_r2 * 5.0;
            }

            else if (TDatabase::ParamDB->stddev_switch == 1)
            {
                double E = TDatabase::ParamDB->stddev_denom;
                double disp = TDatabase::ParamDB->stddev_disp;
                double power = TDatabase::ParamDB->stddev_power;
                double sig_r1 = exp(-pow((2 * actual_x - 1 - disp), power) / (E)) / (2 * 3.14159265359 * sqrt(E)) * exp(-pow((2 * actual_y - 1 - disp), power) / (E)) / (2 * 3.14159265359 * sqrt(E));
                double sig_r2 = exp(-pow((2 * local_x - 1 - disp), power) / (E)) / (2 * 3.14159265359 * sqrt(E)) * exp(-pow((2 * local_y - 1 - disp), power) / (E)) / (2 * 3.14159265359 * sqrt(E));
                // Co Variance
                C[j * N_DOF + i] *= sig_r1 * sig_r2;
            }

            else
            {
                cout << "Error - No standard deviation function is defined for stddev_switch: " << TDatabase::ParamDB->stddev_switch << endl;
                exit(0);
            }

            norm += C[j * N + i] * C[j * N + i];
        }
    }

    ////////////////////////////////////////////////////// SVD ////////////////////////////////////////////
    // Declare SVD parameters
    MKL_INT m1 = N_DOF, n = N_DOF, lda = N_DOF, ldu = N_DOF, ldvt = N_DOF, info;
    double superb[std::min(N_DOF, N_DOF) - 1];

    double *S = new double[N_DOF];
    double *U = new double[N_DOF * N_DOF];
    double *Vt = new double[N_DOF * N_DOF];
    info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', m1, n, C, lda,
                          S, U, ldu, Vt, ldvt, superb);

    cout << endl
         << endl;

    if (info > 0)
    {
        printf("The algorithm computing SVD of Covariance Matrix failed to converge.\n");
        exit(1);
    }
    int energyVal = 0;
    int temp = 0;

    double sumSingularVal = 0;
    for (int i = 0; i < N_DOF; i++)
        sumSingularVal += S[i];

    double val = 0;
    for (energyVal = 0; energyVal < N_DOF; energyVal++)
    {
        val += S[energyVal];
        temp++;
        if (val / sumSingularVal > 0.99)
            break;
    }

    cout << " MODES : " << temp + 1 << endl;

    int modDim = temp + 1;

    double *Ut = new double[N_DOF * modDim]();
    double *Z = new double[N_Realisations * modDim]();

    double *RealizationVector = new double[N_DOF * N_Realisations]();
    double *RealizationVectorTemp = new double[N_DOF * N_Realisations]();

    // -------------- Generate Random Number Based on Normal Distribution -------------------------//
    int k = 0;
    int skip = N_DOF - modDim;
    int count = 0;
    for (int i = 0; i < N_DOF * N_DOF; i++)
    {
        if (count < modDim)
        {
            Ut[k] = U[i];

            count++;
            k++;
        }
        else
        {
            i += skip;
            count = 0;
            i--;
        }
    }

    for (int k = 0; k < modDim; k++)
    {
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<> d{0, 1};

        double *norm1 = new double[N_Realisations];

        for (int n = 0; n < N_Realisations; ++n)
        {
            Z[k * N_Realisations + n] = S[k] * d(gen);
        }
    }

    cout << " N_Realisations : " << N_Realisations << endl;
    cout << " MULT START " << endl;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_DOF, N_Realisations, modDim, 1.0, Ut, modDim, Z, N_Realisations, 0.0, RealizationVector, N_Realisations);
    cout << " MULT DONE " << endl;

    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = 0; j < N_Realisations; j++)
        {
            RealizationVectorTemp[mappingArray[i] * N_Realisations + j] = RealizationVector[j + N_Realisations * i];
        }
    }

    memcpy(RealizationVector, RealizationVectorTemp, N_DOF * N_Realisations * SizeOfDouble);

    cout << N_Realisations << " REALISATIONS COMPUTED " << endl;

    /////////////////////////////////////// -------- END OF REALISATION DATA SETS ------------ ////////////////////////////////////////////////////////////////

    ////////////////////////////////////// -------- START OF DO INITIALIZATION ------------ ////////////////////////////////////////////////////////////////

    double *MeanVector = new double[N_DOF * 1](); // overline{C}_{dof} = \sum_{i=1}^{N_Realisations}(C^{i}_{dof})/N_Realisations
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
    // cout << " DO SVD COMPUTED " << endl;

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

    /////Projection Matrix///////////
    ////////////////////////////////
    double *ProjectionVector = new double[N_Realisations * minDim]();

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N_Realisations, minDim, N_DOF, 1.0, PerturbationVector, N_Realisations, L, minDim, 0.0, ProjectionVector, minDim);

    /// Initialize Coefficient Matrix - First subDim columns of Projection Matrix ////////////////////
    double *CoeffVector = new double[N_Realisations * subDim]();

    for (int i = 0; i < N_Realisations; i++)
    {
        for (int j = 0; j < subDim; j++)
        {
            CoeffVector[j * N_Realisations + i] = ProjectionVector[i * minDim + j]; // CoeffVector in Col Major form
        }
    }

    ////////////Initialize Mode Vector - First subDim columns of Left Singular Vector//////////////////
    double *ModeVector = new double[N_DOF * subDim]();

    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = 0; j < subDim; j++)
        {
            ModeVector[j * N_DOF + i] = L[i * minDim + j]; // ModeVector in Col Major form
        }
    }

    ////////////////////////////////////////////DO - Initialization Ends//////////////////////////////////////
    ///////================================================================================//////////////////

    //======================================================================
    // construct all finite element functions
    //======================================================================
    sol = new double[N_DOF];
    rhs = new double[N_DOF];
    oldrhs = new double[N_DOF];

    memset(sol, 0, N_DOF * SizeOfDouble);
    memset(rhs, 0, N_DOF * SizeOfDouble);
    memset(oldrhs, 0, N_DOF * SizeOfDouble);

    solMean = new double[N_DOF];
    rhsMean = new double[N_DOF];
    old_rhsMean = new double[N_DOF];

    memset(solMean, 0, N_DOF * SizeOfDouble);
    memset(rhsMean, 0, N_DOF * SizeOfDouble);
    memset(old_rhsMean, 0, N_DOF * SizeOfDouble);

    double *solModeAll, *rhsModeAll, *oldsolModeAll, *oldrhsModeAll;
    solModeAll = new double[N_DOF * subDim]();
    rhsModeAll = new double[N_DOF * subDim]();
    oldsolModeAll = new double[N_DOF]();
    oldrhsModeAll = new double[N_DOF]();

    Scalar_FeFunction_Mean = new TFEFunction2D(Scalar_FeSpace, (char *)"C_Mean", (char *)"sol", solMean, N_DOF);

    TFEFunction2D **Scalar_FeFunction_ModeAll = new TFEFunction2D *[subDim];
    for (int s = 0; s < subDim; s++)
    {
        Scalar_FeFunction_ModeAll[s] = new TFEFunction2D(Scalar_FeSpace, (char *)"C_Mode", (char *)"Mode Solution", solModeAll + s * N_DOF, N_DOF);
        Scalar_FeFunction_ModeAll[s]->Interpolate(InitialCondition);
    }

    Scalar_FeFunction_Mean->Interpolate(InitialCondition);
    for (int i = 0; i < N_DOF; i++)
    {
        solMean[i] = MeanVector[i];
    }

    for (int j = 0; j < subDim; j++)
    {
        for (int i = 0; i < N_DOF; i++)
        {
            solModeAll[j * N_DOF + i] = ModeVector[j * N_DOF + i];
        }
    }

    //======================================================================
    // /DO - SystemMatrix construction and solution
    //======================================================================
    // Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
    // Solver: AMG_SOLVE (or) GMG  (or) DIRECT
    TSystemTCD2D *SystemMatrix_Mean = new TSystemTCD2D(Scalar_FeSpace, GALERKIN, DIRECT);
    TSystemTCD2D *SystemMatrix_Mode = new TSystemTCD2D(Scalar_FeSpace, GALERKIN, DIRECT);

    // initilize the system matrix with the functions defined in Example file
    SystemMatrix_Mean->Init(DO_Mean_Equation_Coefficients, BoundCondition, BoundValue);
    SystemMatrix_Mode->Init(DO_Mode_Equation_Coefficients, BoundCondition, BoundValue);

    TSystemTCD2D **SystemMatrix_ModeAll = new TSystemTCD2D *[subDim];
    for (int s = 0; s < subDim; s++)
    {
        SystemMatrix_ModeAll[s] = new TSystemTCD2D(Scalar_FeSpace, GALERKIN, DIRECT);
        SystemMatrix_ModeAll[s]->Init(DO_Mode_Equation_Coefficients, BoundCondition, BoundValue);
    }

    //

    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //------------------------------------------ MEAN EQUATION SETUP -----------------------------------------------------//
    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//

    // Aux Setup for the RHS -- There is no Aux for the Mean equation, So set the values as NULL
    fesp[0] = Scalar_FeSpace;
    aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

    /* -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-
    --------------------------------------[[[ END  ]]] MEAN EQUATION SETUP -----------------------------------------------------
     -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-*/
    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    //  -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //------------------------------------------ CO EFFICIENT EQUATION SETUP -----------------------------------------------------//
    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//

    TFEVectFunct2D *FeVector_Coefficient = new TFEVectFunct2D(Scalar_FeSpace, (char *)"sol", (char *)"sol", CoeffVector, N_Realisations, subDim);

    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //------------------------------------------ CO EFFICIENT EQUATION SETUP END -------------------------------------------------------//
    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//

    // Aux Setup for the RHS -- There is no Aux for the Mean equation, So set the values as NULL
    fesp[0] = Scalar_FeSpace;
    aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //-------------------------------------- MODE EQUATION SETUP -----------------------------------------------------//
    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//

    double *solMode = new double[N_DOF * subDim]();
    double *rhsMode = new double[N_DOF * subDim]();
    double *old_rhsMode = new double[N_DOF]();

    for (int j = 0; j < subDim; j++)
    {
        for (int i = 0; i < N_DOF; i++)
        {
            solMode[j * N_DOF + i] = ModeVector[j * N_DOF + i];
        }
    }

    TFEVectFunct2D *FEFVector_Mode = new TFEVectFunct2D(Scalar_FeSpace, (char *)"C_Mode", (char *)"sol", solMode, N_DOF, subDim);

    int TimeLinear_FESpaces_DO = 1;
    int TimeLinear_Fct_DO = 1; // \tilde(C)
    int TimeLinear_ParamFct_DO = 1;
    int TimeLinear_FEValues_DO = 3;
    int TimeLinear_Params_DO = 3;
    int TimeNSFEFctIndex_DO[3] = {0, 0, 0};
    MultiIndex2D TimeNSFEMultiIndex_DO[3] = {D00, D01, D10};
    ParamFct *TimeNSFct_DO[1] = {DO_Mode_RHS_Aux_Param};
    int TimeNSBeginParam_DO[1] = {0};

    TFEFunction2D *fefct_RHS[4];
    TFESpace2D *fesp_RHS[2];

    TAuxParam2D *aux_RHS_DO = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //--------------------------------------[[[ END  ]]] MODE EQUATION SETUP -----------------------------------------------------//
    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    // assemble the system matrix with given aux, sol and rhs
    // aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
    // otherwise, just pass with NULL
    SystemMatrix_Mean->AssembleMRhs(NULL, solMean, rhsMean);

    SystemMatrix_Mode->AssembleMRhs(NULL, solMode, rhsMode);

    for (int s = 0; s < subDim; s++)
    {
        SystemMatrix_ModeAll[s]->AssembleMRhs(NULL, solModeAll + s * N_DOF, rhsModeAll + s * N_DOF);
    }

    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    //======================================================================
    // produce outout at t=0
    //======================================================================
    // VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
    std::string strMean = std::to_string(N_Realisations);
    std::string filenameMean = "Mean_NRealisations_" + std::to_string(N_Realisations);
    VtkBaseNameMean = const_cast<char *>(filenameMean.c_str());

    std::string fileoutCoeff;
    std::string fileoutMean;
    std::string fileoutMode;
    std::string fileoutMC;

    OutputMean = new TOutput2D(2, 2, 1, 1, Domain);
    OutputMean->AddFEFunction(Scalar_FeFunction_Mean);

    OutputMode = new TOutput2D(2, 2, 1, 1, Domain);
    OutputMode->AddFEVectFunct(FEFVector_Mode);

    TOutput2D **OutputModeAll = new TOutput2D *[subDim];
    for (int s = 0; s < subDim; s++)
    {
        OutputModeAll[s] = new TOutput2D(2, 2, 1, 1, Domain);
        OutputModeAll[s]->AddFEFunction(Scalar_FeFunction_ModeAll[s]);
    }

    int meanimg = 0;
    int modeimg = 0;
    int *imgMode = new int[subDim]();

    if (TDatabase::ParamDB->WRITE_VTK)
    {
        os.seekp(std::ios::beg);
        if (meanimg < 10)
            os << "VTK/" << VtkBaseNameMean << ".0000" << meanimg << ".vtk" << ends;
        else if (meanimg < 100)
            os << "VTK/" << VtkBaseNameMean << ".000" << meanimg << ".vtk" << ends;
        else if (meanimg < 1000)
            os << "VTK/" << VtkBaseNameMean << ".00" << meanimg << ".vtk" << ends;
        else if (meanimg < 10000)
            os << "VTK/" << VtkBaseNameMean << ".0" << meanimg << ".vtk" << ends;
        else
            os << "VTK/" << VtkBaseNameMean << "." << meanimg << ".vtk" << ends;
        OutputMean->WriteVtk(os.str().c_str());
        meanimg++;
    }
    for (int s = 0; s < subDim; s++)
    {
        std::string filenameMode = "Mode_" + std::to_string(s) + "_NRealisations_" + std::to_string(N_Realisations);
        VtkBaseNameMode = const_cast<char *>(filenameMode.c_str());
        if (TDatabase::ParamDB->WRITE_VTK)
        {
            os.seekp(std::ios::beg);
            if (imgMode[s] < 10)
                os << "VTK/" << VtkBaseNameMode << ".0000" << imgMode[s] << ".vtk" << ends;
            else if (imgMode[s] < 100)
                os << "VTK/" << VtkBaseNameMode << ".000" << imgMode[s] << ".vtk" << ends;
            else if (imgMode[s] < 1000)
                os << "VTK/" << VtkBaseNameMode << ".00" << imgMode[s] << ".vtk" << ends;
            else if (imgMode[s] < 10000)
                os << "VTK/" << VtkBaseNameMode << ".0" << imgMode[s] << ".vtk" << ends;
            else
                os << "VTK/" << VtkBaseNameMode << "." << imgMode[s] << ".vtk" << ends;
            OutputModeAll[s]->WriteVtk(os.str().c_str());
            imgMode[s]++;
        }
    }

    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    coll->GetHminHmax(&hmin, &hmax);
    OutPut("h_min : " << hmin << " h_max : " << hmax << endl);

    //======================================================================
    // time disc loop
    //======================================================================
    // parameters for time stepping scheme
    m = 0;
    N_SubSteps = GetN_SubSteps();
    end_time = TDatabase::TimeDB->ENDTIME;

    UpdateStiffnessMat = FALSE; // check BilinearCoeffs in example file
    UpdateRhs = FALSE;          // check BilinearCoeffs in example file
    ConvectionFirstTime = TRUE;
    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    //======================================================================
    // time disc loop
    //======================================================================
    // parameters for time stepping scheme
    m = 0;

    if (m < 10)
        fileoutMean = "Mean/Mean_NRealisations_" + std::to_string(N_Realisations) + "_t0000" + std::to_string(m) + ".txt";
    else if (m < 100)
        fileoutMean = "Mean/Mean_NRealisations_" + std::to_string(N_Realisations) + "_t000" + std::to_string(m) + ".txt";
    else if (m < 1000)
        fileoutMean = "Mean/Mean_NRealisations_" + std::to_string(N_Realisations) + "_t00" + std::to_string(m) + ".txt";
    else if (m < 10000)
        fileoutMean = "Mean/Mean_NRealisations_" + std::to_string(N_Realisations) + "_t0" + std::to_string(m) + ".txt";
    else
        fileoutMean = "Mean/Mean_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m) + ".txt";

    std::ofstream fileMean;
    fileMean.open(fileoutMean);

    for (int i = 0; i < N_DOF; i++)
    {

        fileMean << MeanVector[i] << ",";
        fileMean << endl;
    }

    fileMean.close();

    if (m < 10)
        fileoutMode = "Modes/Mode_NRealisations_" + std::to_string(N_Realisations) + "_t0000" + std::to_string(m) + ".txt";
    else if (m < 100)
        fileoutMode = "Modes/Mode_NRealisations_" + std::to_string(N_Realisations) + "_t000" + std::to_string(m) + ".txt";
    else if (m < 1000)
        fileoutMode = "Modes/Mode_NRealisations_" + std::to_string(N_Realisations) + "_t00" + std::to_string(m) + ".txt";
    else if (m < 10000)
        fileoutMode = "Modes/Mode_NRealisations_" + std::to_string(N_Realisations) + "_t0" + std::to_string(m) + ".txt";
    else
        fileoutMode = "Modes/Mode_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m) + ".txt";
    std::ofstream fileMode;
    fileMode.open(fileoutMode);

    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = 0; j < subDim; j++)
        {
            fileMode << solMode[i * subDim + j];
            if (j != subDim - 1)
                fileMode << ",";
        }
        fileMode << endl;
    }

    fileMode.close();

    if (m < 10)
        fileoutCoeff = "Coefficients/Coeff_NRealisations_" + std::to_string(N_Realisations) + "_t0000" + std::to_string(m) + ".txt";
    else if (m < 100)
        fileoutCoeff = "Coefficients/Coeff_NRealisations_" + std::to_string(N_Realisations) + "_t000" + std::to_string(m) + ".txt";
    else if (m < 1000)
        fileoutCoeff = "Coefficients/Coeff_NRealisations_" + std::to_string(N_Realisations) + "_t00" + std::to_string(m) + ".txt";
    else if (m < 10000)
        fileoutCoeff = "Coefficients/Coeff_NRealisations_" + std::to_string(N_Realisations) + "_t0" + std::to_string(m) + ".txt";
    else
        fileoutCoeff = "Coefficients/Coeff_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m) + ".txt";

    std::ofstream fileCoeff;
    fileCoeff.open(fileoutCoeff);

    for (int i = 0; i < N_Realisations; i++)
    {
        for (int j = 0; j < subDim; j++)
        {
            fileCoeff << CoeffVector[i + j * N_Realisations];
            if (j != subDim - 1)
                fileCoeff << ",";
        }
        fileCoeff << endl;
    }

    fileCoeff.close();

    TDatabase::ParamDB->N_Subspace_Dim = subDim;

    TDatabase::ParamDB->COVARIANCE_MATRIX_DO = new double[subDim * subDim]();

    CalcCovarianceMatx(CoeffVector);

    std::string fileoutMFE;
    std::string fileoutOrtho;
    std::string fileoutprincVar;
    std::ofstream fileMFE;
    std::ofstream fileOrtho;
    std::ofstream fileprincVar;
    double *mfe = new double[(subDim * subDim) + 1]();
    double *princVariances = new double[subDim]();

    memset(mfe, 0, ((subDim * subDim) + 1) * SizeOfDouble);
    calc_MeanFieldEnergy(Scalar_FeSpace, Scalar_FeFunction_Mean, FEFVector_Mode, mfe, subDim);
    

    calc_princVariance(princVariances, subDim);

    if (m < 10)
        fileoutMFE = "Energy_Data/MFE_NRealisations_" + std::to_string(N_Realisations) + "_t0000" + std::to_string(m) + ".txt";
    else if (m < 100)
        fileoutMFE = "Energy_Data/MFE_NRealisations_" + std::to_string(N_Realisations) + "_t000" + std::to_string(m) + ".txt";
    else if (m < 1000)
        fileoutMFE = "Energy_Data/MFE_NRealisations_" + std::to_string(N_Realisations) + "_t00" + std::to_string(m) + ".txt";
    else if (m < 10000)
        fileoutMFE = "Energy_Data/MFE_NRealisations_" + std::to_string(N_Realisations) + "_t0" + std::to_string(m) + ".txt";
    else
        fileoutMFE = "Energy_Data/MFE_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m) + ".txt";

    fileMFE.open(fileoutMFE);

    fileMFE << mfe[0] << endl;
    fileMFE.close();

    if (m < 10)
        fileoutOrtho = "Energy_Data/Ortho_NRealisations_" + std::to_string(N_Realisations) + "_t0000" + std::to_string(m) + ".txt";
    else if (m < 100)
        fileoutOrtho = "Energy_Data/Ortho_NRealisations_" + std::to_string(N_Realisations) + "_t000" + std::to_string(m) + ".txt";
    else if (m < 1000)
        fileoutOrtho = "Energy_Data/Ortho_NRealisations_" + std::to_string(N_Realisations) + "_t00" + std::to_string(m) + ".txt";
    else if (m < 10000)
        fileoutOrtho = "Energy_Data/Ortho_NRealisations_" + std::to_string(N_Realisations) + "_t0" + std::to_string(m) + ".txt";
    else
        fileoutOrtho = "Energy_Data/Ortho_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m) + ".txt";

    fileOrtho.open(fileoutOrtho);
    for (int i = 0; i < subDim; i++)
    {
        for (int j = 0; j < subDim; j++)
        {
            fileOrtho << mfe[i + subDim * j + 1];
            if (j != subDim - 1)
                fileOrtho << ",";
        }
        fileOrtho << endl;
    }

    fileOrtho.close();

    if (m < 10)
        fileoutprincVar = "Energy_Data/PrincVariances_NRealisations_" + std::to_string(N_Realisations) + "_t0000" + std::to_string(m) + ".txt";
    else if (m < 100)
        fileoutprincVar = "Energy_Data/PrincVariances_NRealisations_" + std::to_string(N_Realisations) + "_t000" + std::to_string(m) + ".txt";
    else if (m < 1000)
        fileoutprincVar = "Energy_Data/PrincVariances_NRealisations_" + std::to_string(N_Realisations) + "_t00" + std::to_string(m) + ".txt";
    else if (m < 10000)
        fileoutprincVar = "Energy_Data/PrincVariances_NRealisations_" + std::to_string(N_Realisations) + "_t0" + std::to_string(m) + ".txt";
    else
        fileoutprincVar = "Energy_Data/PrincVariances_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m) + ".txt";

    fileprincVar.open(fileoutprincVar);
    for (int i = 0; i < subDim; i++)
    {
        fileprincVar << princVariances[i] << endl;
    }
    fileprincVar.close();
    // xxxxxx
    N_SubSteps = GetN_SubSteps();
    end_time = TDatabase::TimeDB->ENDTIME;

    UpdateStiffnessMat = TRUE; // check BilinearCoeffs in example file
    UpdateRhs = TRUE;
    ConvectionFirstTime = TRUE;

    // time loop starts
    while (TDatabase::TimeDB->CURRENTTIME < end_time)
    {
        m++;
        TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

        for (l = 0; l < N_SubSteps; l++) // sub steps of fractional step theta
        {
            SetTimeDiscParameters(1);

            if (m == 1)
            {
                OutPut("Theta1: " << TDatabase::TimeDB->THETA1 << endl);
                OutPut("Theta2: " << TDatabase::TimeDB->THETA2 << endl);
                OutPut("Theta3: " << TDatabase::TimeDB->THETA3 << endl);
                OutPut("Theta4: " << TDatabase::TimeDB->THETA4 << endl);
            }

            tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
            TDatabase::TimeDB->CURRENTTIME += tau;

            OutPut(endl
                   << "CURRENT TIME: ");
            OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

            // copy rhs to oldrhs
            memcpy(old_rhsMean, rhsMean, N_DOF * SizeOfDouble);

            // unless the stiffness matrix or rhs change in time, it is enough to
            // assemble only once at the begning
            SystemMatrix_Mean->AssembleARhs(NULL, solMean, rhsMean);

            for (int subSpaceNum = 0; subSpaceNum < subDim; subSpaceNum++)
            {
                // solve the system matrix
                DO_CoEfficient(Scalar_FeSpace, FEFVector_Mode, FeVector_Coefficient, subDim, subSpaceNum, N_Realisations);

                memcpy(old_rhsMode, rhsModeAll + subSpaceNum * N_DOF, N_DOF * SizeOfDouble);

                SystemMatrix_ModeAll[subSpaceNum]->AssembleARhs(NULL, solModeAll + subSpaceNum * N_DOF, rhsModeAll + subSpaceNum * N_DOF);

                // get the UPDATED RHS VALUE FROM FUNCTION
                DO_Mode_RHS(Scalar_FeSpace, FEFVector_Mode, subDim, rhsModeAll + subSpaceNum * N_DOF, subSpaceNum);

                SystemMatrix_ModeAll[subSpaceNum]->AssembleSystMat(old_rhsMode, solModeAll + subSpaceNum * N_DOF, rhsModeAll + subSpaceNum * N_DOF, solModeAll + subSpaceNum * N_DOF);
                SystemMatrix_ModeAll[subSpaceNum]->Solve(solModeAll + subSpaceNum * N_DOF, rhsModeAll + subSpaceNum * N_DOF);
                SystemMatrix_ModeAll[subSpaceNum]->RestoreMassMat();

            } // subSpaceNumLoop

            SystemMatrix_Mean->AssembleSystMat(old_rhsMean, solMean, rhsMean, solMean);
            ////  --
            SystemMatrix_Mean->Solve(solMean, rhsMean);
            SystemMatrix_Mean->RestoreMassMat();

            // restore the mass matrix for the next time step
            // unless the stiffness matrix or rhs change in time, it is not necessary to assemble the system matrix in every time step

        } // for(l=0;l< N_SubSteps;l++)

        //======================================================================
        // produce outout
        //======================================================================

        if (m < 10)
            fileoutMean = "Mean/Mean_NRealisations_" + std::to_string(N_Realisations) + "_t0000" + std::to_string(m) + ".txt";
        else if (m < 100)
            fileoutMean = "Mean/Mean_NRealisations_" + std::to_string(N_Realisations) + "_t000" + std::to_string(m) + ".txt";
        else if (m < 1000)
            fileoutMean = "Mean/Mean_NRealisations_" + std::to_string(N_Realisations) + "_t00" + std::to_string(m) + ".txt";
        else if (m < 10000)
            fileoutMean = "Mean/Mean_NRealisations_" + std::to_string(N_Realisations) + "_t0" + std::to_string(m) + ".txt";
        else
            fileoutMean = "Mean/Mean_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m) + ".txt";

        std::ofstream fileMean;
        fileMean.open(fileoutMean);

        for (int i = 0; i < N_DOF; i++)
        {

            fileMean << MeanVector[i] << ",";
            fileMean << endl;
        }

        fileMean.close();

        if (m < 10)
            fileoutMode = "Modes/Mode_NRealisations_" + std::to_string(N_Realisations) + "_t0000" + std::to_string(m) + ".txt";
        else if (m < 100)
            fileoutMode = "Modes/Mode_NRealisations_" + std::to_string(N_Realisations) + "_t000" + std::to_string(m) + ".txt";
        else if (m < 1000)
            fileoutMode = "Modes/Mode_NRealisations_" + std::to_string(N_Realisations) + "_t00" + std::to_string(m) + ".txt";
        else if (m < 10000)
            fileoutMode = "Modes/Mode_NRealisations_" + std::to_string(N_Realisations) + "_t0" + std::to_string(m) + ".txt";
        else
            fileoutMode = "Modes/Mode_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m) + ".txt";
        std::ofstream fileMode;
        fileMode.open(fileoutMode);

        for (int i = 0; i < N_DOF; i++)
        {
            for (int j = 0; j < subDim; j++)
            {
                fileMode << solMode[i * subDim + j];
                if (j != subDim - 1)
                    fileMode << ",";
            }
            fileMode << endl;
        }

        fileMode.close();

        if (m < 10)
            fileoutCoeff = "Coefficients/Coeff_NRealisations_" + std::to_string(N_Realisations) + "_t0000" + std::to_string(m) + ".txt";
        else if (m < 100)
            fileoutCoeff = "Coefficients/Coeff_NRealisations_" + std::to_string(N_Realisations) + "_t000" + std::to_string(m) + ".txt";
        else if (m < 1000)
            fileoutCoeff = "Coefficients/Coeff_NRealisations_" + std::to_string(N_Realisations) + "_t00" + std::to_string(m) + ".txt";
        else if (m < 10000)
            fileoutCoeff = "Coefficients/Coeff_NRealisations_" + std::to_string(N_Realisations) + "_t0" + std::to_string(m) + ".txt";
        else
            fileoutCoeff = "Coefficients/Coeff_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m) + ".txt";

        std::ofstream fileCoeff;
        fileCoeff.open(fileoutCoeff);

        for (int i = 0; i < N_Realisations; i++)
        {
            for (int j = 0; j < subDim; j++)
            {
                fileCoeff << CoeffVector[i + j * N_Realisations];
                if (j != subDim - 1)
                    fileCoeff << ",";
            }
            fileCoeff << endl;
        }

        fileCoeff.close();

        memset(mfe, 0, ((subDim * subDim) + 1) * SizeOfDouble);
        calc_MeanFieldEnergy(Scalar_FeSpace, Scalar_FeFunction_Mean, FEFVector_Mode, mfe, subDim);
        CalcCovarianceMatx(CoeffVector);

        calc_princVariance(princVariances, subDim);

        if (m < 10)
            fileoutMFE = "Energy_Data/MFE_NRealisations_" + std::to_string(N_Realisations) + "_t0000" + std::to_string(m) + ".txt";
        else if (m < 100)
            fileoutMFE = "Energy_Data/MFE_NRealisations_" + std::to_string(N_Realisations) + "_t000" + std::to_string(m) + ".txt";
        else if (m < 1000)
            fileoutMFE = "Energy_Data/MFE_NRealisations_" + std::to_string(N_Realisations) + "_t00" + std::to_string(m) + ".txt";
        else if (m < 10000)
            fileoutMFE = "Energy_Data/MFE_NRealisations_" + std::to_string(N_Realisations) + "_t0" + std::to_string(m) + ".txt";
        else
            fileoutMFE = "Energy_Data/MFE_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m) + ".txt";

        fileMFE.open(fileoutMFE);

        fileMFE << mfe[0] << endl;
        fileMFE.close();

        if (m < 10)
            fileoutOrtho = "Energy_Data/Ortho_NRealisations_" + std::to_string(N_Realisations) + "_t0000" + std::to_string(m) + ".txt";
        else if (m < 100)
            fileoutOrtho = "Energy_Data/Ortho_NRealisations_" + std::to_string(N_Realisations) + "_t000" + std::to_string(m) + ".txt";
        else if (m < 1000)
            fileoutOrtho = "Energy_Data/Ortho_NRealisations_" + std::to_string(N_Realisations) + "_t00" + std::to_string(m) + ".txt";
        else if (m < 10000)
            fileoutOrtho = "Energy_Data/Ortho_NRealisations_" + std::to_string(N_Realisations) + "_t0" + std::to_string(m) + ".txt";
        else
            fileoutOrtho = "Energy_Data/Ortho_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m) + ".txt";

        fileOrtho.open(fileoutOrtho);
        for (int i = 0; i < subDim; i++)
        {
            for (int j = 0; j < subDim; j++)
            {
                fileOrtho << mfe[i + subDim * j + 1];
                if (j != subDim - 1)
                    fileOrtho << ",";
            }
            fileOrtho << endl;
        }

        fileOrtho.close();

        if (m < 10)
            fileoutprincVar = "Energy_Data/PrincVariances_NRealisations_" + std::to_string(N_Realisations) + "_t0000" + std::to_string(m) + ".txt";
        else if (m < 100)
            fileoutprincVar = "Energy_Data/PrincVariances_NRealisations_" + std::to_string(N_Realisations) + "_t000" + std::to_string(m) + ".txt";
        else if (m < 1000)
            fileoutprincVar = "Energy_Data/PrincVariances_NRealisations_" + std::to_string(N_Realisations) + "_t00" + std::to_string(m) + ".txt";
        else if (m < 10000)
            fileoutprincVar = "Energy_Data/PrincVariances_NRealisations_" + std::to_string(N_Realisations) + "_t0" + std::to_string(m) + ".txt";
        else
            fileoutprincVar = "Energy_Data/PrincVariances_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m) + ".txt";

        fileprincVar.open(fileoutprincVar);
        for (int i = 0; i < subDim; i++)
        {
            fileprincVar << princVariances[i] << endl;
        }
        fileprincVar.close();

        if (m == 1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
        {
            if (TDatabase::ParamDB->WRITE_VTK)
            {
                os.seekp(std::ios::beg);
                if (meanimg < 10)
                    os << "VTK/" << VtkBaseNameMean << ".0000" << meanimg << ".vtk" << ends;
                else if (meanimg < 100)
                    os << "VTK/" << VtkBaseNameMean << ".000" << meanimg << ".vtk" << ends;
                else if (meanimg < 1000)
                    os << "VTK/" << VtkBaseNameMean << ".00" << meanimg << ".vtk" << ends;
                else if (meanimg < 10000)
                    os << "VTK/" << VtkBaseNameMean << ".0" << meanimg << ".vtk" << ends;
                else
                    os << "VTK/" << VtkBaseNameMean << "." << meanimg << ".vtk" << ends;
                OutputMean->WriteVtk(os.str().c_str());
                meanimg++;
            }
            for (int s = 0; s < subDim; s++)
            {
                std::string filenameMode = "Mode_" + std::to_string(s) + "_NRealisations_" + std::to_string(N_Realisations);
                VtkBaseNameMode = const_cast<char *>(filenameMode.c_str());
                if (TDatabase::ParamDB->WRITE_VTK)
                {
                    os.seekp(std::ios::beg);
                    if (imgMode[s] < 10)
                        os << "VTK/" << VtkBaseNameMode << ".0000" << imgMode[s] << ".vtk" << ends;
                    else if (imgMode[s] < 100)
                        os << "VTK/" << VtkBaseNameMode << ".000" << imgMode[s] << ".vtk" << ends;
                    else if (imgMode[s] < 1000)
                        os << "VTK/" << VtkBaseNameMode << ".00" << imgMode[s] << ".vtk" << ends;
                    else if (imgMode[s] < 10000)
                        os << "VTK/" << VtkBaseNameMode << ".0" << imgMode[s] << ".vtk" << ends;
                    else
                        os << "VTK/" << VtkBaseNameMode << "." << imgMode[s] << ".vtk" << ends;
                    OutputModeAll[s]->WriteVtk(os.str().c_str());
                    imgMode[s]++;
                }
            }
        } // produce output loop end
          //======================================================================
        // measure errors to known solution
        //======================================================================
        if (TDatabase::ParamDB->MEASURE_ERRORS)
        {
            fesp[0] = Scalar_FeSpace;
            aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
            Scalar_FeFunction_Mean->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, DO_Mean_Equation_Coefficients, aux, 1, fesp, errors);

            OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
            OutPut(" L2: " << errors[0]);
            OutPut(" H1-semi: " << errors[1] << endl);

            errors[3] += (errors[0] * errors[0] + olderror * olderror) * TDatabase::TimeDB->TIMESTEPLENGTH / 2.0;
            olderror = errors[0];
            OutPut(TDatabase::TimeDB->CURRENTTIME << " L2(0,T;L2) " << sqrt(errors[3]) << " ");

            errors[4] += (errors[1] * errors[1] + olderror1 * olderror1) * TDatabase::TimeDB->TIMESTEPLENGTH / 2.0;
            OutPut("L2(0,T;H1) " << sqrt(errors[4]) << endl);
            olderror1 = errors[1];

            if (Linfty < errors[0])
                Linfty = errors[0];

            OutPut("Linfty " << Linfty << endl);
        } //  if(TDatabase::ParamDB->MEASURE_ERRORS)

    } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

    //======================================================================
    // produce final outout
    //======================================================================
    if (m < 10)
        fileoutMean = "Mean/Mean_NRealisations_" + std::to_string(N_Realisations) + "_t0000" + std::to_string(m) + ".txt";
    else if (m < 100)
        fileoutMean = "Mean/Mean_NRealisations_" + std::to_string(N_Realisations) + "_t000" + std::to_string(m) + ".txt";
    else if (m < 1000)
        fileoutMean = "Mean/Mean_NRealisations_" + std::to_string(N_Realisations) + "_t00" + std::to_string(m) + ".txt";
    else if (m < 10000)
        fileoutMean = "Mean/Mean_NRealisations_" + std::to_string(N_Realisations) + "_t0" + std::to_string(m) + ".txt";
    else
        fileoutMean = "Mean/Mean_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m) + ".txt";

    // std::ofstream fileMean;
    fileMean.open(fileoutMean);

    for (int i = 0; i < N_DOF; i++)
    {

        fileMean << MeanVector[i] << ",";
        fileMean << endl;
    }

    fileMean.close();

    if (m < 10)
        fileoutMode = "Modes/Mode_NRealisations_" + std::to_string(N_Realisations) + "_t0000" + std::to_string(m) + ".txt";
    else if (m < 100)
        fileoutMode = "Modes/Mode_NRealisations_" + std::to_string(N_Realisations) + "_t000" + std::to_string(m) + ".txt";
    else if (m < 1000)
        fileoutMode = "Modes/Mode_NRealisations_" + std::to_string(N_Realisations) + "_t00" + std::to_string(m) + ".txt";
    else if (m < 10000)
        fileoutMode = "Modes/Mode_NRealisations_" + std::to_string(N_Realisations) + "_t0" + std::to_string(m) + ".txt";
    else
        fileoutMode = "Modes/Mode_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m) + ".txt";
    // std::ofstream fileMode;
    fileMode.open(fileoutMode);

    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = 0; j < subDim; j++)
        {
            fileMode << solMode[i * subDim + j];
            if (j != subDim - 1)
                fileMode << ",";
        }
        fileMode << endl;
    }

    fileMode.close();

    if (m < 10)
        fileoutCoeff = "Coefficients/Coeff_NRealisations_" + std::to_string(N_Realisations) + "_t0000" + std::to_string(m) + ".txt";
    else if (m < 100)
        fileoutCoeff = "Coefficients/Coeff_NRealisations_" + std::to_string(N_Realisations) + "_t000" + std::to_string(m) + ".txt";
    else if (m < 1000)
        fileoutCoeff = "Coefficients/Coeff_NRealisations_" + std::to_string(N_Realisations) + "_t00" + std::to_string(m) + ".txt";
    else if (m < 10000)
        fileoutCoeff = "Coefficients/Coeff_NRealisations_" + std::to_string(N_Realisations) + "_t0" + std::to_string(m) + ".txt";
    else
        fileoutCoeff = "Coefficients/Coeff_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m) + ".txt";

    // std::ofstream fileCoeff;
    fileCoeff.open(fileoutCoeff);

    for (int i = 0; i < N_Realisations; i++)
    {
        for (int j = 0; j < subDim; j++)
        {
            fileCoeff << CoeffVector[i + j * N_Realisations];
            if (j != subDim - 1)
                fileCoeff << ",";
        }
        fileCoeff << endl;
    }

    fileCoeff.close();

    memset(mfe, 0, ((subDim * subDim) + 1) * SizeOfDouble);
    calc_MeanFieldEnergy(Scalar_FeSpace, Scalar_FeFunction_Mean, FEFVector_Mode, mfe, subDim);

    CalcCovarianceMatx(CoeffVector);

    calc_princVariance(princVariances, subDim);
    if (m < 10)
        fileoutMFE = "Energy_Data/MFE_NRealisations_" + std::to_string(N_Realisations) + "_t0000" + std::to_string(m) + ".txt";
    else if (m < 100)
        fileoutMFE = "Energy_Data/MFE_NRealisations_" + std::to_string(N_Realisations) + "_t000" + std::to_string(m) + ".txt";
    else if (m < 1000)
        fileoutMFE = "Energy_Data/MFE_NRealisations_" + std::to_string(N_Realisations) + "_t00" + std::to_string(m) + ".txt";
    else if (m < 10000)
        fileoutMFE = "Energy_Data/MFE_NRealisations_" + std::to_string(N_Realisations) + "_t0" + std::to_string(m) + ".txt";
    else
        fileoutMFE = "Energy_Data/MFE_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m) + ".txt";

    fileMFE.open(fileoutMFE);

    fileMFE << mfe[0] << endl;
    fileMFE.close();

    if (m < 10)
        fileoutOrtho = "Energy_Data/Ortho_NRealisations_" + std::to_string(N_Realisations) + "_t0000" + std::to_string(m) + ".txt";
    else if (m < 100)
        fileoutOrtho = "Energy_Data/Ortho_NRealisations_" + std::to_string(N_Realisations) + "_t000" + std::to_string(m) + ".txt";
    else if (m < 1000)
        fileoutOrtho = "Energy_Data/Ortho_NRealisations_" + std::to_string(N_Realisations) + "_t00" + std::to_string(m) + ".txt";
    else if (m < 10000)
        fileoutOrtho = "Energy_Data/Ortho_NRealisations_" + std::to_string(N_Realisations) + "_t0" + std::to_string(m) + ".txt";
    else
        fileoutOrtho = "Energy_Data/Ortho_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m) + ".txt";

    fileOrtho.open(fileoutOrtho);
    for (int i = 0; i < subDim; i++)
    {
        for (int j = 0; j < subDim; j++)
        {
            fileOrtho << mfe[i + subDim * j + 1];
            if (j != subDim - 1)
                fileOrtho << ",";
        }
        fileOrtho << endl;
    }

    fileOrtho.close();
    if (m < 10)
        fileoutprincVar = "Energy_Data/PrincVariances_NRealisations_" + std::to_string(N_Realisations) + "_t0000" + std::to_string(m) + ".txt";
    else if (m < 100)
        fileoutprincVar = "Energy_Data/PrincVariances_NRealisations_" + std::to_string(N_Realisations) + "_t000" + std::to_string(m) + ".txt";
    else if (m < 1000)
        fileoutprincVar = "Energy_Data/PrincVariances_NRealisations_" + std::to_string(N_Realisations) + "_t00" + std::to_string(m) + ".txt";
    else if (m < 10000)
        fileoutprincVar = "Energy_Data/PrincVariances_NRealisations_" + std::to_string(N_Realisations) + "_t0" + std::to_string(m) + ".txt";
    else
        fileoutprincVar = "Energy_Data/PrincVariances_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m) + ".txt";

    fileprincVar.open(fileoutprincVar);
    for (int i = 0; i < subDim; i++)
    {
        fileprincVar << princVariances[i] << endl;
    }
    fileprincVar.close();

    if (TDatabase::ParamDB->WRITE_VTK)
    {
        os.seekp(std::ios::beg);
        if (meanimg < 10)
            os << "VTK/" << VtkBaseNameMean << ".0000" << meanimg << ".vtk" << ends;
        else if (meanimg < 100)
            os << "VTK/" << VtkBaseNameMean << ".000" << meanimg << ".vtk" << ends;
        else if (meanimg < 1000)
            os << "VTK/" << VtkBaseNameMean << ".00" << meanimg << ".vtk" << ends;
        else if (meanimg < 10000)
            os << "VTK/" << VtkBaseNameMean << ".0" << meanimg << ".vtk" << ends;
        else
            os << "VTK/" << VtkBaseNameMean << "." << meanimg << ".vtk" << ends;
        OutputMean->WriteVtk(os.str().c_str());
        meanimg++;
    }
    for (int s = 0; s < subDim; s++)
    {
        std::string filenameMode = "Mode_" + std::to_string(s) + "_NRealisations_" + std::to_string(N_Realisations);
        VtkBaseNameMode = const_cast<char *>(filenameMode.c_str());
        if (TDatabase::ParamDB->WRITE_VTK)
        {
            os.seekp(std::ios::beg);
            if (imgMode[s] < 10)
                os << "VTK/" << VtkBaseNameMode << ".0000" << imgMode[s] << ".vtk" << ends;
            else if (imgMode[s] < 100)
                os << "VTK/" << VtkBaseNameMode << ".000" << imgMode[s] << ".vtk" << ends;
            else if (imgMode[s] < 1000)
                os << "VTK/" << VtkBaseNameMode << ".00" << imgMode[s] << ".vtk" << ends;
            else if (imgMode[s] < 10000)
                os << "VTK/" << VtkBaseNameMode << ".0" << imgMode[s] << ".vtk" << ends;
            else
                os << "VTK/" << VtkBaseNameMode << "." << imgMode[s] << ".vtk" << ends;
            OutputModeAll[s]->WriteVtk(os.str().c_str());
            imgMode[s]++;
        }
    }

    cout << "Subspace Dimension = " << subDim << endl;

    TDatabase::TimeDB->CURRENTTIME = 0;
    std::string PyInFile = "PyIn.txt";

    std::ofstream PyIn;
    PyIn.open(PyInFile);
    PyIn << N_Realisations << endl
         << subDim << endl
         << m << endl;

    CloseFiles();
    return 0;
} // end main
