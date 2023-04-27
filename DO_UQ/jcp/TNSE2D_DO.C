/**
 * @file TNSE2D_ParMooN_IVUQ_DO_Validation.C
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
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <FESpace2D.h>
#include <SystemTNSE2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <Output2D.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>
// #include <TimeUtilities.h>
#include <TNSE2D_ParamRout.h>
#include <TimeDiscRout.h>

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

// =======================================================================
// include utility functions
// =======================================================================
#include "../do_utils/IO.h"                    // Input-Output functions for DO
#include "../do_utils/Monte_Carlo/RealznGen.h" // Monte Carlo Routines
#include "../do_utils/DO_UQ/CoeffOps.h"        //Operations on Coefficient Vector for DO
#include "../do_utils/DO_UQ/ModeOps.h"         //Operations on Mode Vector for DO
// #include "../do_utils/DO_UQ/Stats.h"           //Operations on Mode Vector for DO
#include "../do_utils/DO_UQ/DOInit.h" //Operations on Mode Vector for DO
#include "../do_utils/DO_UQ/DOValidation.h"
#include "../Scratch/DO_Functions.h"
// =======================================================================
// include current example
// =======================================================================

#include "../Examples/TNSE_2D/DrivenCavity.h" //   in unit square

// =======================================================================

int main(int argc, char *argv[])
{

    // ======================================================================
    //  declaration of variables
    // ======================================================================
    int i, j, l, m, N_Cells, ORDER, N_U, N_P, N_L, N_TotalDOF, img = 1, pressure_space_code;
    int Max_It, NSEType, velocity_space_code, N_SubSteps, Disctype;

    double limit, AllErrors[7], end_time, oldtau, tau;

    TDomain *Domain;
    TDatabase *Database = new TDatabase();
    TFEDatabase2D *FEDatabase = new TFEDatabase2D();
    TCollection *coll, *mortarcoll = NULL;
    TFESpace2D *Velocity_FeSpace, *Pressure_FeSpace, *fesp[2];
    TFEVectFunct2D *Velocity;
    TFEFunction2D *u1, *u2, *Pressure, *fefct[2];
    TOutput2D *Output;
    TSystemTNSE2D *SystemMatrix;
    TAuxParam2D *aux, *NSEaux_error;
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
#ifdef __PRIVATE__
    TFESpace2D *Projection_space;
#endif

    const char vtkdir[] = "VTK";
    mkdir(vtkdir, 0777);

    char *PsBaseName, *VtkBaseName, *GEO;
    char UString[] = "u";
    char PString[] = "p";
    char NameString[] = "VMS";

    std::ostringstream os;
    os << " ";

    // ======================================================================
    // set the database values and generate mesh
    // ======================================================================
    /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
    Domain = new TDomain(argv[1]);

    OpenFiles();
    OutFile.setf(std::ios::scientific);

    Database->CheckParameterConsistencyNSE();
    Database->WriteParamDB(argv[0]);
    Database->WriteTimeDB();
    ExampleFile();

    /* include the mesh from a meshgenerator, for a standard mesh use the build-in function */
    // standard mesh
    GEO = TDatabase::ParamDB->GEOFILE;
    Domain->Init(NULL, GEO);

    // refine grid up to the coarsest level
    for (i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
        Domain->RegRefineAll();

    if (TDatabase::ParamDB->WRITE_PS)
    {
        // write grid into an Postscript file
        os.seekp(std::ios::beg);
        os << "Domain"
           << ".ps" << ends;
        Domain->PS(os.str().c_str(), It_Finest, 0);
    }

    //=========================================================================
    // construct all finite element spaces
    //=========================================================================
    ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
    NSEType = TDatabase::ParamDB->NSTYPE;
    Disctype = TDatabase::ParamDB->DISCTYPE;

    coll = Domain->GetCollection(It_Finest, 0);
    N_Cells = coll->GetN_Cells();
    OutPut("N_Cells : " << N_Cells << endl);

    // fespaces for velocity and pressure
    GetVelocityAndPressureSpace(coll, BoundCondition, mortarcoll, Velocity_FeSpace,
                                Pressure_FeSpace, &pressure_space_code,
                                TDatabase::ParamDB->VELOCITY_SPACE,
                                TDatabase::ParamDB->PRESSURE_SPACE);

    // defaulty inf-sup pressure space will be selected based on the velocity space, so update it in database
    TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
    velocity_space_code = TDatabase::ParamDB->VELOCITY_SPACE;

    N_U = Velocity_FeSpace->GetN_DegreesOfFreedom();
    N_P = Pressure_FeSpace->GetN_DegreesOfFreedom();
    N_TotalDOF = 2 * N_U + N_P;

    OutPut("Dof Velocity : " << setw(10) << 2 * N_U << endl);
    OutPut("Dof Pressure : " << setw(10) << N_P << endl);
    OutPut("Total Dof all: " << setw(10) << N_TotalDOF << endl);

#ifdef __PRIVATE__
    if (Disctype == VMS_PROJECTION)
    {
        if (TDatabase::ParamDB->VMS_LARGE_VELOCITY_SPACE == 0)
            Projection_space = new TFESpace2D(coll, NameString, UString, BoundCondition, DiscP_PSpace, 0, mortarcoll);
        else
            Projection_space = new TFESpace2D(coll, NameString, UString, BoundCondition, DiscP_PSpace, 1, mortarcoll);

        N_L = Projection_space->GetN_DegreesOfFreedom();
        OutPut("Dof Projection : " << setw(10) << N_L << endl);
    }
#endif

    ///////////////////////////////////////////////////////////////////////////////////////////////
    ////////// -------- REALISATION DATA GENERATION ----------------------------------------- //////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    int N_Realisations = TDatabase::ParamDB->REALIZATIONS;
    double *RealizationVector = new double[N_U * N_Realisations]();

    GenerateRealizations(Velocity_FeSpace, Pressure_FeSpace, RealizationVector);
    /////////////////////////////////////// -------- END OF REALISATION DATA SETS ------------ ////////////////////////////////////////////////////////////////

    ////////////////////////////////////// -------- START OF DO INITIALIZATION ------------ ////////////////////////////////////////////////////////////////

    double *MeanVector = new double[N_U * 1]();

    int subDim = calculateStochSubspaceDim(Velocity_FeSpace, RealizationVector);

    double *CoeffVector = new double[N_Realisations * subDim]();
    double *ModeVector = new double[N_U * subDim]();
    InitializeDO(Velocity_FeSpace, RealizationVector, MeanVector, ModeVector, CoeffVector);

    /////////////////////////////////////// -------- END OF DO INITIALIZATION ------------ ////////////////////////////////////////////////////////////////

    TFESpace2D *VelocityMean_FeSpace, *PressureMean_FeSpace;
    // fespaces for velocity and pressure
    GetVelocityAndPressureSpace(coll, BoundCondition, mortarcoll, VelocityMean_FeSpace,
                                PressureMean_FeSpace, &pressure_space_code,
                                TDatabase::ParamDB->VELOCITY_SPACE,
                                TDatabase::ParamDB->PRESSURE_SPACE);

    // defaulty inf-sup pressure space will be selected based on the velocity space, so update it in database
    TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
    velocity_space_code = TDatabase::ParamDB->VELOCITY_SPACE;

    int N_U_Mean = VelocityMean_FeSpace->GetN_DegreesOfFreedom();
    int N_P_Mean = PressureMean_FeSpace->GetN_DegreesOfFreedom();
    int N_Total_MeanDOF = 2 * N_U_Mean + N_P_Mean;

    TFESpace2D *VelocityMode_FeSpace, *PressureMode_FeSpace;
    // fespaces for velocity and pressure
    GetVelocityAndPressureSpace(coll, BoundCondition, mortarcoll, VelocityMode_FeSpace,
                                PressureMode_FeSpace, &pressure_space_code,
                                TDatabase::ParamDB->VELOCITY_SPACE,
                                TDatabase::ParamDB->PRESSURE_SPACE);

    // defaulty inf-sup pressure space will be selected based on the velocity space, so update it in database
    TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
    velocity_space_code = TDatabase::ParamDB->VELOCITY_SPACE;

    int N_U_Mode = VelocityMode_FeSpace->GetN_DegreesOfFreedom();
    int N_P_Mode = PressureMode_FeSpace->GetN_DegreesOfFreedom();
    int N_Total_ModeDOF = (2 * N_U_Mode + N_P_Mode) * subDim;

    N_TotalDOF = N_Total_MeanDOF + N_Total_ModeDOF;
    OutPut("Total Mean DOF : 2 * " << setw(10) << N_U_Mean << setw(10) << " + " << setw(10) << N_P_Mean << setw(10) << " = " << setw(10) << N_Total_MeanDOF << endl);
    OutPut("Total Mode DOF : (2 * " << setw(10) << N_U_Mode << setw(10) << " + " << setw(10) << N_P_Mode << setw(10) << ") * " << setw(10) << subDim << setw(10) << " = " << setw(10) << N_Total_ModeDOF << endl);
    OutPut("Total DOF Mean+Mode : " << setw(10) << N_TotalDOF << endl);

    double *solMean, *rhsMean, *oldrhsMean, *solMode, *rhsMode, *oldrhsMode, *sol, *rhs, *oldrhs, *defect, t1, t2, residual, impuls_residual;

    solMean = new double[N_Total_MeanDOF];
    rhsMean = new double[N_Total_MeanDOF];
    oldrhsMean = new double[N_Total_MeanDOF];

    solMode = new double[N_Total_ModeDOF];
    rhsMode = new double[N_Total_ModeDOF];
    oldrhsMode = new double[N_Total_ModeDOF];

    //=========================================================================
    // Assign dimension values to Database
    //=========================================================================
    TDatabase::ParamDB->N_Subspace_Dim = subDim; // Added to Database.h
    TDatabase::ParamDB->REALIZATIONS = N_Realisations;

    //======================================================================
    // ******************** Mean Equation ********************
    //======================================================================

    //======================================================================
    // construct all finite element functions needed for the mean solution
    //======================================================================
    double *solMean, *rhsMean, *oldrhsMean; //,*defect, t1, t2, residual, impuls_residual;
    solMean = new double[N_Total_MeanDOF]();
    rhsMean = new double[N_Total_MeanDOF]();
    oldrhsMean = new double[N_Total_MeanDOF]();

    TFEVectFunct2D *VelocityMean = new TFEVectFunct2D(VelocityMean_FeSpace, (char *)"VelocityMean", (char *)"VelocityMean", solMean, N_U_Mean, 2);
    TFEFunction2D *PressureMean = new TFEFunction2D(PressureMean_FeSpace, (char *)"PressureMean", (char *)"PressureMean", solMean + 2 * N_U_Mean, N_P_Mean);
    TFEFunction2D *u1Mean, *u2Mean;
    u1Mean = VelocityMean->GetComponent(0);
    u2Mean = VelocityMean->GetComponent(1);

    for (i = 0; i < N_U_Mean; i++)
    {

        solMean[i] = MeanVector[i];
        solMean[N_U_Mean + i] = MeanVector[i];
    }

    for (i = 0; i < N_P_Mean; i++)
        solMean[2 * N_U_Mean + i] = 0;

    //======================================================================
    // ******************** Mode Equation ********************
    //======================================================================

    //======================================================================
    // construct all finite element functions needed for the mode solution
    //======================================================================
    double *solMode, *rhsMode, *old_rhsMode, *defectMode, *old_solMode;
    solMode = new double[N_Total_ModeDOF]();
    old_solMode = new double[2 * N_U_Mode + N_P_Mode]();
    rhsMode = new double[N_Total_ModeDOF]();
    old_rhsMode = new double[2 * N_U_Mode + N_P_Mode]();

    for (int j = 0; j < subDim; j++)
    {
        for (int i = 0; i < N_U_Mode; i++)
        {
            solMode[(j * (2 * N_U_Mode + N_P_Mode)) + i] = ModeVector[j * N_U_Mode + i];            // first component velocity
            solMode[(j * (2 * N_U_Mode + N_P_Mode)) + N_U_Mode + i] = ModeVector[j * N_U_Mode + i]; // second component velocity
        }
    }

    for (int j = 0; j < subDim; j++)
    {
        for (int i = 0; i < N_P_Mode; i++)
        {

            solMode[(j * (2 * N_U_Mode + N_P_Mode)) + N_U_Mode + N_U_Mode + i] = 0; // Pressure
        }
    }

    double *solModeVeloCopy = new double[2 * N_U_Mode * subDim]();

    for (int i = 0; i < subDim; i++)
    {
        for (int j = 0; j < 2 * N_U_Mode; j++)
            solModeVeloCopy[2 * N_U_Mode * i + j] = solMode[i * (2 * N_U_Mode + N_P_Mode) + j];
    }

    double *solModePressCopy = new double[N_P_Mode * subDim]();

    for (int i = 0; i < subDim; i++)
    {
        for (int j = 0; j < N_P_Mode; j++)
            solModePressCopy[N_P_Mode * i + j] = solMode[i * (2 * N_U_Mode + N_P_Mode) + 2 * N_U_Mode + j];
        0
    }

    TFEVectFunct2D *VelocityMode = new TFEVectFunct2D(VelocityMode_FeSpace, (char *)"VelocityMode", (char *)"VelocityMode", solModeVeloCopy, N_U_Mode, 2 * subDim);

    TFEFunction2D *PressureMode = new TFEFunction2D(PressureMode_FeSpace, (char *)"PressureMode", (char *)"PressureMode", solModePressCopy, N_P_Mode);

    TFEFunction2D *u1Mode, *u2Mode;
    u1Mode = VelocityMode->GetComponent(0);
    u2Mode = VelocityMode->GetComponent(1);

    TFEVectFunct2D **VelocityModeAll = new TFEVectFunct2D *[subDim];
    for (int subD = 0; subD < subDim; subD++)
        VelocityModeAll[subD] = new TFEVectFunct2D(VelocityMode_FeSpace, (char *)"U_Mode", (char *)"Mode Component", solMode + (subD * (2 * N_U_Mode + N_P_Mode)), N_U_Mode, 2);

    TFEFunction2D **PressureModeAll = new TFEFunction2D *[subDim];
    for (int subD = 0; subD < subDim; subD++)
        PressureModeAll[subD] = new TFEFunction2D(PressureMode_FeSpace, (char *)"P_Mode", (char *)"Mode Component", solMode + (subDim * (2 * N_U_Mode + N_P_Mode)) + 2 * N_U_Mode, N_P_Mode);

    //======================================================================
    // SystemMatrix construction for Mean Equation
    //======================================================================

    TSystemTNSE2D *SystemMatrixMean = new TSystemTNSE2D(VelocityMean_FeSpace, PressureMean_FeSpace, VelocityMean, PressureMean, solMean, rhsMean, Disctype, NSEType, DIRECT
#ifdef __PRIVATE__
                                                        ,
                                                        Projection_space, NULL, NULL
#endif
    );

    TFESpace2D *fespMean[2];

    fespMean[0] = VelocityMean_FeSpace;

    TFEFunction2D *fefctMean[2];
    fefctMean[0] = u1Mean;
    fefctMean[1] = u2Mean;

    TAuxParam2D *auxMean, *NSEaux_error_mean;

    switch (Disctype)
    {
    // turbulent viscosity must be computed
    case SMAGORINSKY:
    case VMS_PROJECTION:
    case CLASSICAL_LES:
    case GL00_CONVOLUTION:
    case GL00_AUX_PROBLEM:

        auxMean = new TAuxParam2D(TimeNSN_FESpacesVelo_GradVelo, TimeNSN_FctVelo_GradVelo,
                                  TimeNSN_ParamFctVelo_GradVelo,
                                  TimeNSN_FEValuesVelo_GradVelo,
                                  fespMean, fefctMean,
                                  TimeNSFctVelo_GradVelo,
                                  TimeNSFEFctIndexVelo_GradVelo,
                                  TimeNSFEMultiIndexVelo_GradVelo,
                                  TimeNSN_ParamsVelo_GradVelo,
                                  TimeNSBeginParamVelo_GradVelo);

        break;

    default:
        // 2 parameters are needed for assembling (u1_old, u2_old)
        auxMean = new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2, TimeNSN_ParamFct2,
                                  TimeNSN_FEValues2,
                                  fespMean, fefctMean,
                                  TimeNSFct2,
                                  TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
                                  TimeNSN_Params2, TimeNSBeginParam2);
    }

    // aux for calculating the error
    if (TDatabase::ParamDB->MEASURE_ERRORS)
    {
        NSEaux_error_mean = new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
                                            TimeNSN_ParamFct2,
                                            TimeNSN_FEValues2,
                                            fespMean, fefctMean,
                                            TimeNSFct2,
                                            TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
                                            TimeNSN_Params2, TimeNSBeginParam2);
    }

    TAuxParam2D *auxMode, *NSEaux_error_mode;
    TFESpace2D *fespMode[2];
    TFEFunction2D *fefctModeAll[50][2];

    fespMode[0] = VelocityMode_FeSpace;

    TAuxParam2D **auxModeAll = new TAuxParam2D *[subDim];
    TAuxParam2D **NSEaux_error_modeAll = new TAuxParam2D *[subDim];

    for (int subD = 0; subD < subDim; subD++)
    {
        fefctModeAll[subD][0] = VelocityModeAll[subD]->GetComponent(0);
        fefctModeAll[subD][1] = VelocityModeAll[subD]->GetComponent(1);
        switch (Disctype)
        {
        // turbulent viscosity must be computed
        case SMAGORINSKY:
        case VMS_PROJECTION:
        case CLASSICAL_LES:
        case GL00_CONVOLUTION:
        case GL00_AUX_PROBLEM:

            auxModeAll[subD] = new TAuxParam2D(TimeNSN_FESpacesVelo_GradVelo, TimeNSN_FctVelo_GradVelo,
                                               TimeNSN_ParamFctVelo_GradVelo,
                                               TimeNSN_FEValuesVelo_GradVelo,
                                               fespMode, fefctModeAll[subD],
                                               TimeNSFctVelo_GradVelo,
                                               TimeNSFEFctIndexVelo_GradVelo,
                                               TimeNSFEMultiIndexVelo_GradVelo,
                                               TimeNSN_ParamsVelo_GradVelo,
                                               TimeNSBeginParamVelo_GradVelo);

            break;

        default:
            // 2 parameters are needed for assembling (u1_old, u2_old)
            auxModeAll[subD] = new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2, TimeNSN_ParamFct2,
                                               TimeNSN_FEValues2,
                                               fespMode, fefctModeAll[subD],
                                               TimeNSFct2,
                                               TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
                                               TimeNSN_Params2, TimeNSBeginParam2);
        }

        // aux for calculating the error
        if (TDatabase::ParamDB->MEASURE_ERRORS)
        {
            NSEaux_error_modeAll[subD] = new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
                                                         TimeNSN_ParamFct2,
                                                         TimeNSN_FEValues2,
                                                         fespMode, fefctModeAll[subD],
                                                         TimeNSFct2,
                                                         TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
                                                         TimeNSN_Params2, TimeNSBeginParam2);
        }
    }

    return 0;
}