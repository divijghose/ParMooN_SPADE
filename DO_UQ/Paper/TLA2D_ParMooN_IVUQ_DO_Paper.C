/**
 * @file TLA2D_ParMooN_IVUQ_DO_Validation.C
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

 *@bug Stochastic Normalization of Modes for inner product -
       The inner product in the finite element calculation doesn't produce the same inner product matrix as the vector sense, the modes have to be normalized for this *Divij - 01/08/22*
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
#include "../HelperFunctions/IO.h"                    // Input-Output functions for DO
#include "../HelperFunctions/Monte_Carlo/RealznGen.h" // Monte Carlo Routines
#include "../HelperFunctions/DO_UQ/CoeffOps.h"        //Operations on Coefficient Vector for DO
#include "../HelperFunctions/DO_UQ/ModeOps.h"         //Operations on Mode Vector for DO
// #include "../HelperFunctions/DO_UQ/Stats.h"           //Operations on Mode Vector for DO
#include "../HelperFunctions/DO_UQ/DOInit.h" //Operations on Mode Vector for DO
#include "../HelperFunctions/DO_UQ/DOValidation.h"

#include "../../Examples/DO_UQ/linear_advection_do_validation.h"

int main(int argc, char *argv[])
{
    int i, j, l, m, N_SubSteps, ORDER, N_Cells, N_DOF, img = 1, N_G;
    int N_Active;

    double *sol, *rhs, *oldrhs, t1, t2, errors[5], Linfty;
    double tau, end_time, *defect, olderror, olderror1, hmin, hmax;

    bool UpdateStiffnessMat, UpdateRhs, ConvectionFirstTime;

    char *VtkBaseName;
    char *VtkBaseNameMean;
    char *VtkBaseNameMode;

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

    // create output directories, if not already existing
    createOutputFolders();

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

    int N_Realisations = TDatabase::ParamDB->REALIZATIONS;
    double *RealizationVector = new double[N_DOF * N_Realisations]();

    GenerateRealizations(Scalar_FeSpace, RealizationVector);

    ////////////////////////////////////// -------- START OF DO INITIALIZATION ------------ ////////////////////////////////////////////////////////////////
    double *MeanVector = new double[N_DOF * 1]();

    int subDim = calculateStochSubspaceDim(Scalar_FeSpace, RealizationVector);

    double *CoeffVector = new double[N_Realisations * subDim]();
    double *ModeVector = new double[N_DOF * subDim]();
    InitializeDO(Scalar_FeSpace, RealizationVector, MeanVector, ModeVector, CoeffVector);

    ////////////////////////////////////// -------- END OF DO INITIALIZATION ------------ ////////////////////////////////////////////////////////////////

    double *IPMatxMode = new double[subDim * subDim]();
    double *IPMatxMean = new double[1 * 1]();
    m = 0;

    //======================================================================
    // construct all finite element functions
    //======================================================================
    sol = new double[N_DOF]();
    rhs = new double[N_DOF]();
    oldrhs = new double[N_DOF]();

    double *solMean = new double[N_DOF]();
    double *rhsMean = new double[N_DOF]();
    double *old_rhsMean = new double[N_DOF]();

    double *solModeAll = new double[N_DOF * subDim]();
    double *rhsModeAll = new double[N_DOF * subDim]();
    double *oldsolModeAll = new double[N_DOF]();
    double *oldrhsModeAll = new double[N_DOF]();

    Scalar_FeFunction_Mean = new TFEFunction2D(Scalar_FeSpace, (char *)"C_Mean", (char *)"Mean Solution", solMean, N_DOF);

    TFEFunction2D **Scalar_FeFunction_ModeAll = new TFEFunction2D *[subDim];
    for (int s = 0; s < subDim; s++)
    {
        Scalar_FeFunction_ModeAll[s] = new TFEFunction2D(Scalar_FeSpace, (char *)"C_Mode", (char *)"Mode Solution", solModeAll + s * N_DOF, N_DOF);
    }

    memcpy(solMean, MeanVector, N_DOF * SizeOfDouble);

    for (int j = 0; j < subDim; j++)
    {
        for (int i = 0; i < N_DOF; i++)
        {
            memcpy(solModeAll + j * N_DOF, ModeVector + j * N_DOF, N_DOF * SizeOfDouble);
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

    // Aux Setup for the RHS -- There is no Aux for the Mean equation, So set the values as NULL
    fesp[0] = Scalar_FeSpace;
    aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

    TFEVectFunct2D *FeVector_Coefficient = new TFEVectFunct2D(Scalar_FeSpace, (char *)"Phi", (char *)"Coefficients", CoeffVector, N_Realisations, subDim);

    // Aux Setup for the RHS -- There is no Aux for the Mean equation, So set the values as NULL
    fesp[0] = Scalar_FeSpace;
    aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    double *solMode = new double[N_DOF * subDim]();
    double *rhsMode = new double[N_DOF * subDim]();
    double *old_rhsMode = new double[N_DOF]();

    for (int j = 0; j < subDim; j++)
    {
        for (int i = 0; i < N_DOF; i++)
        {
            memcpy(solModeAll + j * N_DOF, ModeVector + j * N_DOF, N_DOF * SizeOfDouble);
        }
    }

    TFEVectFunct2D *FEFVector_Mode = new TFEVectFunct2D(Scalar_FeSpace, (char *)"C_Mode", (char *)"Mode Solution", solMode, N_DOF, subDim);

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

    
    int *imgMode = new int[subDim]();
    int *imgMean = new int(0);

    printVTKOutput(VtkBaseNameMean, imgMean, OutputMean);

    for (int s = 0; s < subDim; s++)
    {
        std::string filenameMode = "Mode_" + std::to_string(s) + "_NRealisations_" + std::to_string(N_Realisations);
        VtkBaseNameMode = const_cast<char *>(filenameMode.c_str());
        printVTKOutput(VtkBaseNameMode, imgMode + s, OutputModeAll[s]);
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

    //======================================================================
    // time disc loop
    //======================================================================
    // parameters for time stepping scheme
    m = 0;
    std::string fileoutCoeff;
    std::string fileoutMean;
    std::string fileoutMode;
    std::string meanBaseName = "Mean/Mean_NRealisations_";
    std::string modeBaseName = "Modes/Mode_NRealisations_";
    std::string coeffBaseName = "Coefficients/Coeff_NRealisations_";
    std::string IPModeBaseName = "IPMatrices/Mode/IPMode_NRealisations_";
    std::string IPMeanBaseName = "IPMatrices/Mean/IPMean_NRealisations_";

    fileoutMean = generateFileName(meanBaseName, m, N_Realisations);
    printToTxt(fileoutMean, solMean, N_DOF, 1, 'C');

    fileoutMode = generateFileName(modeBaseName, m, N_Realisations);
    printToTxt(fileoutMode, solModeAll, N_DOF, subDim, 'C');

    fileoutCoeff = generateFileName(coeffBaseName, m, N_Realisations);
    printToTxt(fileoutCoeff, CoeffVector, N_Realisations, subDim, 'C');

    std::string fileoutIPMode;
    std::string fileoutIPMean;

    fileoutIPMode = generateFileName(IPModeBaseName, m, N_Realisations);
    calcIPMatx(IPMatxMode, solModeAll, N_DOF, subDim, 'C');
    printToTxt(fileoutIPMode, IPMatxMode, subDim, subDim, 'R');

    fileoutIPMean = generateFileName(IPMeanBaseName, m, N_Realisations);
    calcIPMatx(IPMatxMean, solMean, N_DOF, 1, 'C');
    printToTxt(fileoutIPMean, IPMatxMean, 1, 1, 'R');

    TDatabase::ParamDB->N_Subspace_Dim = subDim;

    TDatabase::ParamDB->COVARIANCE_MATRIX_DO = new double[subDim * subDim]();

    CalcCovarianceMatx(CoeffVector);

    std::string fileoutMFE;
    std::string fileoutOrtho;
    std::string fileoutprincVar;
    std::ofstream fileMFE;
    std::ofstream fileOrtho;
    std::ofstream fileprincVar;

    std::string orthoBaseName = "Energy_Data/IPModeFE/Ortho_NRealisations_";
    std::string mfeBaseName = "Energy_Data/MFE/MFE_NRealisations_";
    std::string princVarBaseName = "Energy_Data/PV/PrincVariances_NRealisations_";

    double *modeOrtho = new double[subDim * subDim]();
    double *princVariances = new double[subDim]();
    double mfe = 0;

    // memset(modeOrtho, 0, (subDim * subDim) * SizeOfDouble);
    mfe = calc_MeanFieldEnergy(Scalar_FeSpace, Scalar_FeFunction_Mean);

    calc_princVariance(princVariances, subDim);

    fileoutMFE = generateFileName(mfeBaseName, m, N_Realisations);
    printToTxt(fileoutMFE, &mfe, 1, 1, 'R');

    fileoutOrtho = generateFileName(orthoBaseName, m, N_Realisations);
    calc_ModeOrtho(Scalar_FeSpace, FEFVector_Mode, subDim, modeOrtho);
    printToTxt(fileoutOrtho, modeOrtho, subDim, subDim, 'R');

    fileoutprincVar = generateFileName(princVarBaseName, m, N_Realisations);
    printToTxt(fileoutprincVar, princVariances, subDim, 1, 'R');

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
            SystemMatrix_Mean->AssembleSystMat(old_rhsMean, solMean, rhsMean, solMean);
            ////  --
            SystemMatrix_Mean->Solve(solMean, rhsMean);
            SystemMatrix_Mean->RestoreMassMat();

            for (int subSpaceNum = 0; subSpaceNum < subDim; subSpaceNum++)
            {
                // solve the system matrix

                memcpy(old_rhsMode, rhsModeAll + subSpaceNum * N_DOF, N_DOF * SizeOfDouble);
                memcpy(oldsolModeAll, solModeAll + subSpaceNum * N_DOF, N_DOF * SizeOfDouble);

                SystemMatrix_ModeAll[subSpaceNum]->AssembleARhs(NULL, solModeAll + subSpaceNum * N_DOF, rhsModeAll + subSpaceNum * N_DOF);
                double *rhsNewAll = new double[N_DOF * subDim]();
                // get the UPDATED RHS VALUE FROM FUNCTION

                DO_Mode_RHS(Scalar_FeSpace, FEFVector_Mode, subDim, rhsModeAll + subSpaceNum * N_DOF, subSpaceNum);

                SystemMatrix_ModeAll[subSpaceNum]->AssembleSystMat(old_rhsMode, solModeAll + subSpaceNum * N_DOF, rhsModeAll + subSpaceNum * N_DOF, solModeAll + subSpaceNum * N_DOF);

                SystemMatrix_ModeAll[subSpaceNum]->Solve(solModeAll + subSpaceNum * N_DOF, rhsModeAll + subSpaceNum * N_DOF);

                SystemMatrix_ModeAll[subSpaceNum]->RestoreMassMat();

                DO_CoEfficient(Scalar_FeSpace, FEFVector_Mode, FeVector_Coefficient, subDim, subSpaceNum, N_Realisations);

            } // subSpaceNumLoop

            // restore the mass matrix for the next time step
            // unless the stiffness matrix or rhs change in time, it is not necessary to assemble the system matrix in every time step

        } // for(l=0;l< N_SubSteps;l++)

        reorthonormalizeC(solModeAll, N_DOF, subDim);

        //======================================================================
        // produce outout
        //======================================================================

        fileoutMean = generateFileName(meanBaseName, m, N_Realisations);
        printToTxt(fileoutMean, solMean, N_DOF, 1, 'C');

        fileoutMode = generateFileName(modeBaseName, m, N_Realisations);
        printToTxt(fileoutMode, solModeAll, N_DOF, subDim, 'C');

        fileoutCoeff = generateFileName(coeffBaseName, m, N_Realisations);
        printToTxt(fileoutCoeff, CoeffVector, N_Realisations, subDim, 'C');

        memset(modeOrtho, 0, (subDim * subDim) * SizeOfDouble);
        mfe = calc_MeanFieldEnergy(Scalar_FeSpace, Scalar_FeFunction_Mean);
        CalcCovarianceMatx(CoeffVector);

        calc_princVariance(princVariances, subDim);

        fileoutMFE = generateFileName(mfeBaseName, m, N_Realisations);
        printToTxt(fileoutMFE, &mfe, 1, 1, 'R');

        fileoutOrtho = generateFileName(orthoBaseName, m, N_Realisations);
        calc_ModeOrtho(Scalar_FeSpace, FEFVector_Mode, subDim, modeOrtho);
        printToTxt(fileoutOrtho, modeOrtho, subDim, subDim, 'R');

        fileoutprincVar = generateFileName(princVarBaseName, m, N_Realisations);
        printToTxt(fileoutprincVar, princVariances, subDim, 1, 'R');

        fileoutIPMode = generateFileName(IPModeBaseName, m, N_Realisations);
        calcIPMatx(IPMatxMode, solModeAll, N_DOF, subDim, 'C');
        printToTxt(fileoutIPMode, IPMatxMode, subDim, subDim, 'R');

        fileoutIPMean = generateFileName(IPMeanBaseName, m, N_Realisations);
        calcIPMatx(IPMatxMean, solMean, N_DOF, 1, 'C');
        printToTxt(fileoutIPMean, IPMatxMean, 1, 1, 'R');

        if (m == 1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
        {
            printVTKOutput(VtkBaseNameMean, imgMean, OutputMean);

            for (int s = 0; s < subDim; s++)
            {
                std::string filenameMode = "Mode_" + std::to_string(s) + "_NRealisations_" + std::to_string(N_Realisations);
                VtkBaseNameMode = const_cast<char *>(filenameMode.c_str());
                printVTKOutput(VtkBaseNameMode, imgMode + s, OutputModeAll[s]);
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

        memcpy(solMode, solModeAll, N_DOF * subDim * SizeOfDouble);
    } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

    //======================================================================
    // produce final outout
    //======================================================================
    m++;
    fileoutMean = generateFileName(meanBaseName, m, N_Realisations);
    printToTxt(fileoutMean, solMean, N_DOF, 1, 'C');

    fileoutMode = generateFileName(modeBaseName, m, N_Realisations);
    printToTxt(fileoutMode, solModeAll, N_DOF, subDim, 'C');

    fileoutCoeff = generateFileName(coeffBaseName, m, N_Realisations);
    printToTxt(fileoutCoeff, CoeffVector, N_Realisations, subDim, 'C');

    memset(modeOrtho, 0, (subDim * subDim) * SizeOfDouble);
    mfe = calc_MeanFieldEnergy(Scalar_FeSpace, Scalar_FeFunction_Mean);

    CalcCovarianceMatx(CoeffVector);

    calc_princVariance(princVariances, subDim);

    fileoutMFE = generateFileName(mfeBaseName, m, N_Realisations);
    printToTxt(fileoutMFE, &mfe, 1, 1, 'R');

    fileoutOrtho = generateFileName(orthoBaseName, m, N_Realisations);
    calc_ModeOrtho(Scalar_FeSpace, FEFVector_Mode, subDim, modeOrtho);
    printToTxt(fileoutOrtho, modeOrtho, subDim, subDim, 'R');

    fileoutprincVar = generateFileName(princVarBaseName, m, N_Realisations);
    printToTxt(fileoutprincVar, princVariances, subDim, 1, 'R');

    fileoutIPMode = generateFileName(IPModeBaseName, m, N_Realisations);
    calcIPMatx(IPMatxMode, solModeAll, N_DOF, subDim, 'C');
    printToTxt(fileoutIPMode, IPMatxMode, subDim, subDim, 'R');

    calcIPMatx(IPMatxMean, solMean, N_DOF, 1, 'C');
    printToTxt(fileoutIPMean, IPMatxMean, 1, 1, 'R');

    printVTKOutput(VtkBaseNameMean, imgMean, OutputMean);

    for (int s = 0; s < subDim; s++)
    {
        std::string filenameMode = "Mode_" + std::to_string(s) + "_NRealisations_" + std::to_string(N_Realisations);
        VtkBaseNameMode = const_cast<char *>(filenameMode.c_str());
        printVTKOutput(VtkBaseNameMode, imgMode + s, OutputModeAll[s]);
    }

    TDatabase::TimeDB->CURRENTTIME = 0;
    std::string PyInFile = "PyIn.txt";

    std::ofstream PyIn;
    PyIn.open(PyInFile);
    PyIn << N_Realisations << endl
         << subDim << endl
         << m << endl;

    CloseFiles();
    TDatabase::TimeDB->CURRENTTIME = 0;

    cout << "Compute finish" << endl;

    Output = new TOutput2D(2, 2, 1, 1, Domain);
    SystemMatrix = new TSystemTCD2D(Scalar_FeSpace, GALERKIN, DIRECT);

    // // initilize the system matrix with the functions defined in Example file
    SystemMatrix->Init(BilinearCoeffs, BoundCondition, BoundValue);
    validateDOvsMC(Scalar_FeSpace, RealizationVector, Output, SystemMatrix);

    TDatabase::TimeDB->CURRENTTIME = 0;

    CloseFiles();
    return 0;

} // end main
