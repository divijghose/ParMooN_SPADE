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
    int N_Realisations = TDatabase::ParamDB->REALISATIONS;
    double *RealizationVector = new double[N_U * N_Realisations]();

    // Generate Realisation Data Sets
    GenerateRealizations(Velocity_FeSpace, Pressure_FeSpace, RealizationVector);

    /////////////////////////////////////// -------- END OF REALISATION DATA SETS ------------ ////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////// Monte-Carlo Simulation ///////////////////////////////////////////////////////////////////////////////////////
    double *solMC, *rhsMC;
    solMC = new double[N_TotalDOF]();
    rhsMC = new double[N_TotalDOF]();


    TFEVectFunct2D *VelocityMC = new TFEVectFunct2D(Velocity_FeSpace, UString, UString, solMC, N_U, 2);
    u1 = VelocityMC->GetComponent(0);
    u2 = VelocityMC->GetComponent(1);

    TFEFunction2D *PressureMC = new TFEFunction2D(Pressure_FeSpace, PString, PString, solMC + (2 * N_U), N_P);

    // define the system matrix
    TSystemTNSE2D *SystemMatrixMC = new TSystemTNSE2D(Velocity_FeSpace, Pressure_FeSpace, VelocityMC, PressureMC, solMC, rhsMC, Disctype, NSEType, DIRECT
#ifdef __PRIVATE__
                                                      ,
                                                      Projection_space, NULL, NULL
#endif
    );

    

    // define the aux
    fesp[0] = Velocity_FeSpace;
    fefct[0] = u1;
    fefct[1] = u2;

    switch (Disctype)
    {
    // turbulent viscosity must be computed
    case SMAGORINSKY:
    case VMS_PROJECTION:
    case CLASSICAL_LES:
    case GL00_CONVOLUTION:
    case GL00_AUX_PROBLEM:

        aux = new TAuxParam2D(TimeNSN_FESpacesVelo_GradVelo, TimeNSN_FctVelo_GradVelo,
                              TimeNSN_ParamFctVelo_GradVelo,
                              TimeNSN_FEValuesVelo_GradVelo,
                              fesp, fefct,
                              TimeNSFctVelo_GradVelo,
                              TimeNSFEFctIndexVelo_GradVelo,
                              TimeNSFEMultiIndexVelo_GradVelo,
                              TimeNSN_ParamsVelo_GradVelo,
                              TimeNSBeginParamVelo_GradVelo);

        break;

    default:
        // 2 parameters are needed for assembling (u1_old, u2_old)
        aux = new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2, TimeNSN_ParamFct2,
                              TimeNSN_FEValues2,
                              fesp, fefct,
                              TimeNSFct2,
                              TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
                              TimeNSN_Params2, TimeNSBeginParam2);
    }

    // aux for calculating the error
    if (TDatabase::ParamDB->MEASURE_ERRORS)
    {
        NSEaux_error = new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
                                       TimeNSN_ParamFct2,
                                       TimeNSN_FEValues2,
                                       fesp, fefct,
                                       TimeNSFct2,
                                       TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
                                       TimeNSN_Params2, TimeNSBeginParam2);
    }
    double hmin, hmax;
    coll->GetHminHmax(&hmin, &hmax);

    // initilize the system matrix with the functions defined in Example file
    // last argument is aux that is used to pass additional fe functions (eg. mesh velocity)
    SystemMatrixMC->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue, aux, NSEaux_error);

    for (int RealNo = 0; RealNo < N_Realisations; RealNo++)
    {
        OutPut("Realisation Number : " << RealNo << endl);
        for (int i = 0; i < N_U; i++)
        {
            solMC[i] = RealizationVector[i * N_Realisations + RealNo];
            // solMC[N_U + i] = RealizationVector[i * N_Realisations + RealNo];
            solMC[N_U + i] = 0;
        }
        SystemMatrixMC->Assemble(solMC, rhsMC);

        mkdir("MC_Realisations", 0777);
        std::string folderNameMC = "MC_Realisations/MC_Realisation_" + std::to_string(RealNo);
        mkdir(folderNameMC.c_str(), 0777);
        
        std::string fileNameMC = "MC_Realisation_" + std::to_string(RealNo);
        char *vtkBaseNameMC = const_cast<char *>(fileNameMC.c_str());


        //======================================================================
        // time disc loop
        //======================================================================
        // parameters for time stepping scheme
        int timeStepCounterMC = 0;
        int numSubSteps = 2;
        double oldTau = 1.;
        double endTimeMC = 0.01;
        limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
        Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
        memset(AllErrors, 0, 7. * SizeOfDouble);

        TOutput2D *OutputMC = new TOutput2D(2, 2, 1, 1, Domain);
        OutputMC->AddFEFunction(PressureMC);
        OutputMC->AddFEVectFunct(VelocityMC);

        int *imgMC = new int(0);

        printVTKOutput(vtkBaseNameMC, imgMC, OutputMC);
    }

    //======================================================================

    ////////////////////////// End of Monte-Carlo Simulation ///////////////////////////////////////////////////////////////////////////////////////

    return 0;
}