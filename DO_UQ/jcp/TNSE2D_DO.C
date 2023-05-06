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
    double *solMC, *rhsMC, *oldrhsMC, *defectMC, residualMC, impuls_residualMC;
    solMC = new double[N_TotalDOF]();
    rhsMC = new double[N_TotalDOF]();
    oldrhsMC = new double[N_TotalDOF]();
    defectMC = new double[N_TotalDOF]();

    TFEVectFunct2D *VelocityMC = new TFEVectFunct2D(Velocity_FeSpace, UString, UString, solMC, N_U, 2);
    u1 = VelocityMC->GetComponent(0);
    u2 = VelocityMC->GetComponent(1);
    //  interpolate the initial solution

    TFEFunction2D *PressureMC = new TFEFunction2D(Pressure_FeSpace, PString, PString, solMC + (2 * N_U), N_P);
    u1->Interpolate(InitialU1);
    u2->Interpolate(InitialU2);
    PressureMC->Interpolate(InitialP);
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

    mkdir("MC_Realisations", 0777);
    std::string folderNameMC;

    std::string fileNameMC;
    char *vtkBaseNameMC;
    //======================================================================
    // time disc loop
    //======================================================================
    // parameters for time stepping scheme
    int timeStepCounterMC = 0;
    int numSubSteps = GetN_SubSteps();
    double oldTau = 1.;
    double endTimeMC = TDatabase::TimeDB->ENDTIME;
    // double endTimeMC = 0.01;
    limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
    Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
    TOutput2D *OutputMC = new TOutput2D(2, 2, 1, 1, Domain);
    OutputMC->AddFEFunction(PressureMC);
    OutputMC->AddFEVectFunct(VelocityMC);

// commenting out 302 to 466

//     for (int RealNo = 0; RealNo < N_Realisations; RealNo++)
//     {
//         // OutPut("Realisation Number : " << RealNo << endl);
//         folderNameMC = "MC_Realisations/MC_Realisation_" + std::to_string(RealNo);
//         mkdir(folderNameMC.c_str(), 0777);

//         for (int i = 0; i < N_U; i++)
//         {
//             solMC[i] = RealizationVector[i * N_Realisations + RealNo];
//             solMC[N_U + i] = RealizationVector[i * N_Realisations + RealNo];
//             // solMC[N_U + i] = 0;
//         }
//         SystemMatrixMC->Assemble(solMC, rhsMC);

//         fileNameMC = "MC_Realisation_" + std::to_string(RealNo);
//         vtkBaseNameMC = const_cast<char *>(fileNameMC.c_str());

//         memset(AllErrors, 0, 7. * SizeOfDouble);

//         int *imgMC = new int(0);

//         printVTKOutput(vtkBaseNameMC, imgMC, OutputMC);

//         // Checked till here
//                // time loop starts while (TDatabase::TimeDB->CURRENTTIME < endTimeMC)
//         { // time cycle
//             timeStepCounterMC++;
//             TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
//             for (l = 0; l < numSubSteps; l++) // sub steps of fractional step theta
//             {
//                 SetTimeDiscParameters(1);

//                 if (timeStepCounterMC == 1)
//                 {
//                     OutPut("Theta1: " << TDatabase::TimeDB->THETA1 << endl);
//                     OutPut("Theta2: " << TDatabase::TimeDB->THETA2 << endl);
//                     OutPut("Theta3: " << TDatabase::TimeDB->THETA3 << endl);
//                     OutPut("Theta4: " << TDatabase::TimeDB->THETA4 << endl);
//                 }

//                 tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
//                 TDatabase::TimeDB->CURRENTTIME += tau;

//                 OutPut(endl
//                        << "CURRENT TIME: ");
//                 OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

//                 // copy sol, rhs to olssol, oldrhs
//                 memcpy(oldrhsMC, rhsMC, N_TotalDOF * SizeOfDouble);

//                 // assemble only rhs, nonlinear matrix for NSE will be assemble in fixed point iteration
//                 // not needed if rhs is not time-dependent
//                 if (timeStepCounterMC != 1)
//                 {
//                     SystemMatrixMC->AssembleRhs(solMC, rhsMC);
//                 }
//                 else
//                 {
//                     SystemMatrixMC->Assemble(solMC, rhsMC);
//                 }

//                 // scale B matices and assemble NSE-rhs based on the \theta time stepping scheme
//                 SystemMatrixMC->AssembleSystMat(tau / oldtau, oldrhsMC, rhsMC, solMC);
//                 oldtau = tau;

//                 // calculate the residual
//                 defectMC = new double[N_TotalDOF];
//                 memset(defectMC, 0, N_TotalDOF * SizeOfDouble);

//                 SystemMatrixMC->GetTNSEResidual(solMC, defectMC);

//                 // correction due to L^2_O Pressure space
//                 if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
//                     IntoL20Vector2D(defectMC + 2 * N_U, N_P, pressure_space_code);

//                 residualMC = Ddot(N_TotalDOF, defectMC, defectMC);
//                 impuls_residualMC = Ddot(2 * N_U, defectMC, defectMC);
//                 OutPut("Nonlinear iteration step   0");
//                 OutPut(setw(14) << impuls_residualMC);
//                 OutPut(setw(14) << residualMC - impuls_residualMC);
//                 OutPut(setw(14) << sqrt(residualMC) << endl);

//                 //======================================================================
//                 // Solve the system
//                 // Nonlinear iteration of fixed point type
//                 //======================================================================
//                 for (j = 1; j <= Max_It; j++)
//                 {
//                     // Solve the NSE system
//                     SystemMatrixMC->Solve(solMC);

//                     if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
//                         IntoL20FEFunction(solMC + 2 * N_U, N_P, Pressure_FeSpace, velocity_space_code, pressure_space_code);

//                     // no nonlinear iteration for Stokes problem
//                     if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == STOKES)
//                         break;

//                     // restore the mass matrix for the next nonlinear iteration
//                     SystemMatrixMC->RestoreMassMat();

//                     // assemble the system matrix with given aux, sol and rhs
//                     SystemMatrixMC->AssembleANonLinear(solMC, rhsMC);

//                     // assemble system mat, S = M + dt\theta_1*A
//                     SystemMatrixMC->AssembleSystMatNonLinear();

//                     // get the residual
//                     memset(defectMC, 0, N_TotalDOF * SizeOfDouble);
//                     SystemMatrixMC->GetTNSEResidual(solMC, defectMC);

//                     // correction due to L^2_O Pressure space
//                     if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
//                         IntoL20Vector2D(defectMC + 2 * N_U, N_P, pressure_space_code);

//                     residualMC = Ddot(N_TotalDOF, defectMC, defectMC);
//                     impuls_residualMC = Ddot(2 * N_U, defectMC, defectMC);
//                     OutPut("nonlinear iteration step " << setw(3) << j);
//                     OutPut(setw(14) << impuls_residualMC);
//                     OutPut(setw(14) << residualMC - impuls_residualMC);
//                     OutPut(setw(14) << sqrt(residualMC) << endl);

//                     if (sqrt(residualMC) <= limit)
//                         break;

//                 } // for(j=1;j<=Max_It;j++)
//                   /*           cout << " test VHM main " << endl;
// exit(0);      */
//                 // restore the mass matrix for the next time step
//                 SystemMatrixMC->RestoreMassMat();

//             } // for(l=0;l<N_SubSteps;
//               //======================================================================
//               // measure errors to known solution
//               //======================================================================
//             if (TDatabase::ParamDB->MEASURE_ERRORS)
//             {
//                 //        u1->Interpolate(ExactU1);
//                 //        u2->Interpolate(ExactU2);
//                 //        Pressure->Interpolate(ExactP);

//                 SystemMatrixMC->MeasureTNSEErrors(ExactU1, ExactU2, ExactP, AllErrors);

//                 OutPut("L2(u): " << AllErrors[0] << endl);
//                 OutPut("H1-semi(u): " << AllErrors[1] << endl);
//                 OutPut("L2(p): " << AllErrors[2] << endl);
//                 OutPut("H1-semi(p): " << AllErrors[3] << endl);
//                 OutPut(AllErrors[4] << " l_infty(L2(u)) " << AllErrors[5] << endl);
//                 OutPut(TDatabase::TimeDB->CURRENTTIME << " L2(0,t,L2)(u) : " << sqrt(AllErrors[6]) << endl);

//             } // if(TDatabase::ParamDB->MEASURE_ERRORS)

//             //======================================================================
//             // produce outout
//             //======================================================================
//             printVTKOutput(vtkBaseNameMC, imgMC, OutputMC);

//         } // while(TDatabase::TimeDB->CURRENTTIME< e
//         delete[] imgMC;
//         TDatabase::TimeDB->CURRENTTIME = 0;
//         //======================================================================
//         // produce final outout
//         //======================================================================

//     } // realisation Loop

    //======================================================================

    ////////////////////////// End of Monte-Carlo Simulation ///////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////Start of DO - Initialization ///////////////////////////////////////////////////////////////////////////////////////
    double *MeanVector = new double[N_U * 1]();

    int subDim = calculateStochSubspaceDim(Velocity_FeSpace, RealizationVector);

    double *CoeffVector = new double[N_Realisations * subDim]();
    double *ModeVector = new double[N_U * subDim]();
    InitializeDO(Velocity_FeSpace, RealizationVector, MeanVector, ModeVector, CoeffVector);
    cout << "Works till Initialize" << endl;

    /////////////////////////////End of DO - Initialization //////////////////////////////////////
    ///////================================================================================//////////////////
    ///////////////////////////////////////////DO - Initialization Ends//////////////////////////////////////
    ///////================================================================================//////////////////

    CalcCovarianceMatx(CoeffVector); // Added to Database.h

    //=========================================================================
    // Assign dimension values to Database
    //=========================================================================
    TDatabase::ParamDB->N_Subspace_Dim = subDim; // Added to Database.h
    TDatabase::ParamDB->REALISATIONS = N_Realisations;

    //=========================================================================
    // Set up data structures for velocity and pressure
    //=========================================================================
    double *MeanVectorComp1 = MeanVector; // Assign Calculated Mean vector to first component of mean velocity (u bar)
    double *MeanVectorComp2 = new double[N_U];
    memset(MeanVectorComp2, 0, N_U * SizeOfDouble); // Second component of mean velocity has been initialized to zero **Has to be changed if uncertainty exists in both components**

    double *ModeVectorComp1 = ModeVector; // Assign calculated Mode vector to mode of first component of velocity (u Tilde)
    double *ModeVectorComp2 = new double[N_U * subDim];
    memset(ModeVectorComp2, 0, N_U * subDim * SizeOfDouble); // Mode of second component of velocity is initialized to zero (v Tilde)** has to be changed if both components have uncertainty**

    //*Aggregate array for velocity mode - Combining the two components of mode vector to form a VectFunction array
    //*Structure of Aggregate array = [[U_Tilde_1 transposed],[V_Tilde_1 transposed],...,...,[U_Tilde_N_S transposed],[V_Tilde_N_S transposed]]
    double *TotalModeVector = new double[N_U * subDim * 2]();




    // Column Major forms of Mode Vector Components
    double *ModeVect1Col = RowtoColMaj(ModeVectorComp1, N_U, subDim);
    double *ModeVect2Col = RowtoColMaj(ModeVectorComp2, N_U, subDim);

    // Total Mode Vector in Column Major Form
    for (int i = 0; i < subDim; i++)
    {
        memcpy(TotalModeVector + (2 * i * N_U), ModeVect1Col + (1 * N_U), N_U * SizeOfDouble);
        memcpy(TotalModeVector + ((2 * i + 1) * N_U), ModeVect2Col + (1 * N_U), N_U * SizeOfDouble);
    }

    //======================================================================
    // construct all finite element functions
    //======================================================================

    // for Mean
    double *solMean, *rhsMean, *old_rhsMean, *defectMean, *old_solMean;
    int N_Total_MeanDOF = 2 * N_U + N_P;
    solMean = new double[N_Total_MeanDOF]();
    rhsMean = new double[N_Total_MeanDOF]();
    old_rhsMean = new double[N_Total_MeanDOF]();
    old_solMean = new double[N_Total_MeanDOF]();


    for (i = 0; i < N_U; i++)
    {

        solMean[i] = MeanVector[i];
        solMean[N_U + i] = MeanVector[i];
        // solMean[i] = 0;
        // solMean[N_U + i] = 0;
        // solMean[N_U + i] = 0;
    }

    for (i = 0; i < N_P; i++)
    {

        solMean[2 * N_U + i] = 0;
    }

    TFEVectFunct2D *Velocity_Mean;
    TFEFunction2D *Pressure_Mean;

    Velocity_Mean = new TFEVectFunct2D(Velocity_FeSpace, (char *)"U_Mean", (char *)"Mean Component", solMean, N_U, 2); // check length

    Pressure_Mean = new TFEFunction2D(Pressure_FeSpace, (char *)"P_Mean", (char *)"Mean Component", solMean + 2 * N_U, N_P);
        // For Mean 
    TFEFunction2D* u1Mean = Velocity_Mean->GetComponent(0);
    TFEFunction2D* u2Mean = Velocity_Mean->GetComponent(1);

    // Interpolate the initial solution
    //  u1Mean->Interpolate(InitialU1);
    //  u2Mean->Interpolate(InitialU2);
    //  Pressure_Mean->Interpolate(InitialP);



    // for mode
    double *solMode, *rhsMode, *old_rhsMode, *defectMode, *old_solMode;
    int N_Total_ModeDOF = (2 * N_U + N_P) * subDim;
    solMode = new double[N_Total_ModeDOF]();
    old_solMode = new double[2 * N_U + N_P]();
    rhsMode = new double[N_Total_ModeDOF]();
    old_rhsMode = new double[2 * N_U + N_P]();

    double *solMode1 = new double[N_U * subDim]();
    double *solMode2 = new double[N_U * subDim]();

    double *solModeAll = new double[N_Total_ModeDOF]();

    for (int j = 0; j < subDim; j++)
    {
        for (int i = 0; i < N_U; i++)
        {
            solMode[(j * (2 * N_U + N_P)) + i] = ModeVector[j * N_U + i]; // first component velocity
            solModeAll[(j * (2 * N_U + N_P)) + i] = ModeVector[j * N_U + i];
            solModeAll[(j * (2 * N_U + N_P)) + i + N_U] = ModeVector[j * N_U + i];

            // solMode[(j * (2 * N_M + N_P)) + N_M + i] = 0;				  // second component velocity
        }
    }

    // for (int j = 0; j < subDim; j++)
    // {
    // 	for (int i = 0; i < N_P; i++)
    // 	{

    // 		solMode[(j * (2 * N_M + N_P)) + N_M + N_M + i] = 0; // Pressure
    // 	}
    // }
    double *solModeVeloCopy = new double[2 * N_U * subDim]();

    for (int i = 0; i < subDim; i++)
    {
        for (int j = 0; j < 2 * N_U; j++)
            solModeVeloCopy[2 * N_U * i + j] = solMode[i * (2 * N_U + N_P) + j];
    }

    double *solModePressCopy = new double[N_P * subDim]();

    for (int i = 0; i < subDim; i++)
    {
        for (int j = 0; j < N_P; j++)
            solModePressCopy[N_P * i + j] = solMode[i * (2 * N_U + N_P) + 2 * N_U + j];
    }

    TFEVectFunct2D *Velocity_Mode;
    TFEFunction2D *u1Mode, *u2Mode, *Pressure_Mode, *fefctMode[2];




    
    ///////////////


    Velocity_Mode = new TFEVectFunct2D(Velocity_FeSpace, (char *)"U_Mode", (char *)"Mode Component", solModeVeloCopy, N_U, 2 * subDim); // check length ??

    // u1Mode = Velocity_Mode->GetComponent(0);
    // u2Mode = Velocity_Mode->GetComponent(1);
    Pressure_Mode = new TFEFunction2D(Pressure_FeSpace, (char *)"P_Mode", (char *)"Mode Component", solModePressCopy, N_P); // check length ??

    // Variables to capture Vorticity for all Modes
    double *vort_mode = new double[N_U * subDim]();
    double *div_mode = new double[N_U * subDim]();

    TFEFunction2D** vorticityFeFunctionMode = new TFEFunction2D*[subDim]();
    for (int i = 0; i < subDim; i++)
    {
        vorticityFeFunctionMode[i] = new TFEFunction2D(Velocity_FeSpace,(char*)"Vorticity",(char*)"Vorticity",vort_mode + i * N_U,N_U);
    }

    TFEVectFunct2D **VelocityModeAll = new TFEVectFunct2D *[subDim];
    for (int subD = 0; subD < subDim; subD++)
    {
        VelocityModeAll[subD] = new TFEVectFunct2D(Velocity_FeSpace, (char *)"U_Mode", (char *)"Mode Component", solModeAll + (subD * (2 * N_U + N_P)), N_U, 2);

    }

    TFEFunction2D **PressureModeAll = new TFEFunction2D *[subDim];
    for (int subD = 0; subD < subDim; subD++)
        PressureModeAll[subD] = new TFEFunction2D(Pressure_FeSpace, (char *)"P_Mode", (char *)"Mode Component", solModeAll + (subDim * (2 * N_U + N_P)) + 2 * N_U, N_P);

    //======================================================================
    // SystemMatrix construction and solution
    //======================================================================

    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //------------------------------------------ MEAN EQUATION SETUP -----------------------------------------------------//
    TSystemTNSE2D *SystemMatrix_Mean;
    SystemMatrix_Mean = new TSystemTNSE2D(Velocity_FeSpace, Pressure_FeSpace, Velocity_Mean, Pressure_Mean, solMean, rhsMean, Disctype, NSEType, DIRECT
#ifdef __PRIVATE__
                                          ,
                                          Projection_space, NULL, NULL
#endif
    );

    // Output Functions
    // Setup Output for Mean and N_S outputs for modes
    char *PsBaseNameMean, *VtkBaseNameMean, *GEOMean;
    std::string strMean = std::to_string(N_Realisations);
    std::string filenameMean = "Mean_NRealisations_" + std::to_string(N_Realisations);
    VtkBaseNameMean = const_cast<char *>(filenameMean.c_str());

    char *PsBaseNameMode, *VtkBaseNameMode, *GEOMode;
    std::string strMode = std::to_string(N_Realisations);
    std::string filenameMode = "Mode_NRealisations_" + std::to_string(N_Realisations);
    VtkBaseNameMode = const_cast<char *>(filenameMode.c_str());

    TOutput2D *OutputMean;
    OutputMean = new TOutput2D(2, 2, 1, 1, Domain);
    OutputMean->AddFEVectFunct(Velocity_Mean);
    OutputMean->AddFEFunction(Pressure_Mean);
    

    TOutput2D **OutputModeAll = new TOutput2D *[subDim];
    for (int subD = 0; subD < subDim; subD++)
    {
        OutputModeAll[subD] = new TOutput2D(2, 2, 1, 1, Domain);
        OutputModeAll[subD]->AddFEVectFunct(VelocityModeAll[subD]);
        OutputModeAll[subD]->AddFEFunction(PressureModeAll[subD]);
    }

    int meanimg = 0;
    int *modeimg = new int[subDim]();

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

    for (int sd = 0; sd < subDim; sd++)
    {

        std::string filenameMode = "Mode_NRealisations_" + std::to_string(N_Realisations) + "ModeN0_" + std::to_string(sd);
        VtkBaseNameMode = const_cast<char *>(filenameMode.c_str());
        if (TDatabase::ParamDB->WRITE_VTK)
        {
            os.seekp(std::ios::beg);
            if (modeimg[sd] < 10)
                os << "VTK/" << VtkBaseNameMode << ".0000" << modeimg[sd] << ".vtk" << ends;
            else if (modeimg[sd] < 100)
                os << "VTK/" << VtkBaseNameMode << ".000" << modeimg[sd] << ".vtk" << ends;
            else if (modeimg[sd] < 1000)
                os << "VTK/" << VtkBaseNameMode << ".00" << modeimg[sd] << ".vtk" << ends;
            else if (modeimg[sd] < 10000)
                os << "VTK/" << VtkBaseNameMode << ".0" << modeimg[sd] << ".vtk" << ends;
            else
                os << "VTK/" << VtkBaseNameMode << "." << modeimg[sd] << ".vtk" << ends;
            OutputModeAll[sd]->WriteVtk(os.str().c_str());
            modeimg[sd]++;
        }
    }
    // cout << "Before First Time Step" << endl;
    // exit(0);
    

    TFESpace2D *fespMean[2];

    fespMean[0] = Velocity_FeSpace;

    TFEFunction2D *fefctMean[2];
    u1Mean = Velocity_Mean->GetComponent(0);
    u2Mean = Velocity_Mean->GetComponent(1);
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
    // double hmin, hmax;
    // coll->GetHminHmax(&hmin, &hmax);
    // OutPut("h_min : " << hmin << " h_max : " << hmax << endl);

    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //------------------------------------------ MEAN EQUATION SETUP END-----------------------------------------------------//

    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //------------------------------------------ MODE EQUATION SETUP -----------------------------------------------------//
    double *rhsModeAll = new double[N_Total_ModeDOF]();
    TSystemTNSE2D **SystemMatrixModeAll = new TSystemTNSE2D *[subDim];
    for (int subD = 0; subD < subDim; subD++)
        SystemMatrixModeAll[subD] = new TSystemTNSE2D(Velocity_FeSpace, Pressure_FeSpace, VelocityModeAll[subD], PressureModeAll[subD], solModeAll + (subD * (2 * N_U + N_P)), rhsModeAll + (subD * (2 * N_U + N_P)), Disctype, NSEType, DIRECT
#ifdef __PRIVATE__
                                                      ,
                                                      Projection_space, NULL, NULL
#endif
        );

    TAuxParam2D *auxMode, *NSEaux_error_mode;
    TFESpace2D *fespMode[2];
    TFEFunction2D *fefctModeAll[50][2];

    fespMode[0] = Velocity_FeSpace;

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

    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //------------------------------------------ MODE EQUATION SETUP END-----------------------------------------------------//
    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //------------------------------------------ COEFF EQUATION SETUP END-----------------------------------------------------//

    TFEVectFunct2D *FeVector_Coefficient = new TFEVectFunct2D(Velocity_FeSpace, (char *)"Coeff", (char *)"Coefficients", CoeffVector, N_Realisations, subDim);
    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //------------------------------------------ COEFF EQUATION SETUP END END-----------------------------------------------------//

    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //------------------------------------------ MEAN EQUATION INITIALIZATION -----------------------------------------------------//
    SystemMatrix_Mean->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue, auxMean, NSEaux_error_mean);
    SystemMatrix_Mean->Assemble(solMean, rhsMean);

    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //------------------------------------------ MODE EQUATION INITIALIZATION -----------------------------------------------------//
    for (int subD = 0; subD < subDim; subD++)
    {
        SystemMatrixModeAll[subD]->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue, auxModeAll[subD], NSEaux_error_modeAll[subD]);
        SystemMatrixModeAll[subD]->Assemble(solModeAll + (subD * (2 * N_U + N_P)), rhsModeAll + (subD * (2 * N_U + N_P))); // seg fault
    }
    TDatabase::ParamDB->COVARIANCE_MATRIX_DO = new double[subDim * subDim]();
    TDatabase::ParamDB->COSKEWNESS_MATRIX_DO = new double[subDim * subDim * subDim]();
    TDatabase::ParamDB->COVARIANCE_INVERSE_DO = new double[subDim * subDim]();

    CalcCovarianceMatx(CoeffVector);
    CalcCoskewnessMatx(CoeffVector);
    InvertCov();

    TDatabase::TimeDB->CURRENTTIME = 0.0;

    // done here -- Commentedby thivin , since duplicate assembly
    // SystemMatrix_Mean->Assemble(solMean, rhsMean);
    // for (int subD = 0; subD < subDim; subD++)
    //     SystemMatrixModeAll[subD]->Assemble(solModeAll + (subD * (2 * N_U + N_P)), rhsModeAll + (subD * (2 * N_U + N_P)));


    // Print after Assembly Routine
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

    for (int sd = 0; sd < subDim; sd++)
    {

        std::string filenameMode = "Mode_NRealisations_" + std::to_string(N_Realisations) + "ModeN0_" + std::to_string(sd);
        VtkBaseNameMode = const_cast<char *>(filenameMode.c_str());
        if (TDatabase::ParamDB->WRITE_VTK)
        {
            os.seekp(std::ios::beg);
            if (modeimg[sd] < 10)
                os << "VTK/" << VtkBaseNameMode << ".0000" << modeimg[sd] << ".vtk" << ends;
            else if (modeimg[sd] < 100)
                os << "VTK/" << VtkBaseNameMode << ".000" << modeimg[sd] << ".vtk" << ends;
            else if (modeimg[sd] < 1000)
                os << "VTK/" << VtkBaseNameMode << ".00" << modeimg[sd] << ".vtk" << ends;
            else if (modeimg[sd] < 10000)
                os << "VTK/" << VtkBaseNameMode << ".0" << modeimg[sd] << ".vtk" << ends;
            else
                os << "VTK/" << VtkBaseNameMode << "." << modeimg[sd] << ".vtk" << ends;
            OutputModeAll[sd]->WriteVtk(os.str().c_str());
            modeimg[sd]++;
        }
    }

    
    TDatabase::TimeDB->CURRENTTIME = 0;
    //======================================================================
    // time disc loop
    //======================================================================
    // parameters for time stepping scheme
    m = 0;
    N_SubSteps = GetN_SubSteps();
    oldtau = 1.;
    end_time = TDatabase::TimeDB->ENDTIME;
    limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
    Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
    memset(AllErrors, 0, 7. * SizeOfDouble);
    defectMean = new double[N_Total_MeanDOF]();
    double residual, impuls_residual;

   
    // time loop starts
    while (TDatabase::TimeDB->CURRENTTIME < end_time)
    { // time cycle
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

            // copy sol, rhs to olssol, oldrhs
            memcpy(old_rhsMean, rhsMean, N_Total_MeanDOF * SizeOfDouble);
            memcpy(old_solMean, solMean, N_Total_MeanDOF * SizeOfDouble);

            // assemble only rhs, nonlinear matrix for NSE will be assemble in fixed point iteration
            // not needed if rhs is not time-dependent
            if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == STOKES)
                DO_Mean_RHS(Velocity_FeSpace, Velocity_Mode, subDim, rhsMean, N_U);

            SystemMatrix_Mean->Assemble(solMean, rhsMean);

            // scale B matices and assemble NSE-rhs based on the \theta time stepping scheme
            SystemMatrix_Mean->AssembleSystMat(tau / oldtau, old_rhsMean, rhsMean, solMean);
            oldtau = tau;

            // calculate the residual

            memset(defectMean, 0, N_Total_MeanDOF * SizeOfDouble);

            SystemMatrix_Mean->GetTNSEResidual(solMean, defectMean);

            // correction due to L^2_O Pressure space
            if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
                IntoL20Vector2D(defectMean + 2 * N_U, N_P, pressure_space_code);

            residual = Ddot(N_Total_MeanDOF, defectMean, defectMean);
            impuls_residual = Ddot(2 * N_U, defectMean, defectMean);
            OutPut("Nonlinear iteration step   0");
            OutPut(setw(14) << impuls_residual);
            OutPut(setw(14) << residual - impuls_residual);
            OutPut(setw(14) << sqrt(residual) << endl);

            //======================================================================
            // Solve the system
            // Nonlinear iteration of fixed point type
            //======================================================================
            for (j = 1; j <= Max_It; j++)
            {
                // Solve the NSE system
                SystemMatrix_Mean->Solve(solMean);

                if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
                    IntoL20FEFunction(solMean + 2 * N_U, N_P, Pressure_FeSpace, velocity_space_code, pressure_space_code);

                // no nonlinear iteration for Stokes problem
                if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == STOKES)
                    break;

                // restore the mass matrix for the next nonlinear iteration
                SystemMatrix_Mean->RestoreMassMat();

                // assemble the system matrix with given aux, sol and rhs
                SystemMatrix_Mean->AssembleANonLinear(solMean, rhsMean);

                // assemble system mat, S = M + dt\theta_1*A
                SystemMatrix_Mean->AssembleSystMatNonLinear();

                // get the residual
                // calculate the residual
                memset(defectMean, 0, N_Total_MeanDOF * SizeOfDouble);

                SystemMatrix_Mean->GetTNSEResidual(solMean, defectMean);

                // correction due to L^2_O Pressure space
                if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
                    IntoL20Vector2D(defectMean + 2 * N_U, N_P, pressure_space_code);

                residual = Ddot(N_Total_MeanDOF, defectMean, defectMean);
                impuls_residual = Ddot(2 * N_U, defectMean, defectMean);
                OutPut("nonlinear iteration step " << setw(3) << j);
                OutPut(setw(14) << impuls_residual);
                OutPut(setw(14) << residual - impuls_residual);
                OutPut(setw(14) << sqrt(residual) << endl);

                if (sqrt(residual) <= limit)
                    break;

            } // for(j=1;j<=Max_It;j++)
              /*           cout << " test VHM main " << endl;
exit(0);      */
            // restore the mass matrix for the next time step
            SystemMatrix_Mean->RestoreMassMat();

            for (int subSpaceNum = 0; subSpaceNum < subDim; subSpaceNum++)
            {

                DO_CoEfficient(Velocity_FeSpace, Velocity_Mode, FeVector_Coefficient, Velocity_Mean, Pressure_Mode, subDim, subSpaceNum, N_Realisations);
            }

            CalcCovarianceMatx(CoeffVector);
            CalcCoskewnessMatx(CoeffVector);
            InvertCov();

            for (int subSpaceNum = 0; subSpaceNum < subDim; subSpaceNum++)
            { // Subspace loop
                // subspace loop
                double *modeSolution_i = solMode + subSpaceNum * (2 * N_U + N_P); // This works for column major
                double *modeSolution_rhs = rhsMode + subSpaceNum * (2 * N_U + N_P);
                // copy sol, rhs to olssol, oldrhs
                memcpy(old_rhsMode, modeSolution_rhs, (2 * N_U + N_P) * SizeOfDouble);
                // memcpy(old_solMode, modeSolution_i, N_Total_MeanDOF * SizeOfDouble);

                // Assemble RHS
                DO_Mode_RHS(Velocity_FeSpace, Velocity_Mean, Velocity_Mode, Pressure_Mode, subDim, modeSolution_rhs, subSpaceNum);

                SystemMatrixModeAll[subSpaceNum]->Assemble(modeSolution_i, modeSolution_rhs);
                // scale B matices and assemble NSE-rhs based on the \theta time stepping scheme

                SystemMatrixModeAll[subSpaceNum]->AssembleSystMat(tau / oldtau, old_rhsMode, modeSolution_rhs, modeSolution_i);
                oldtau = tau;

                // calculate the residual
                defectMode = new double[2 * N_U + N_P]();

                memset(defectMode, 0, (2 * N_U + N_P) * SizeOfDouble);

                SystemMatrixModeAll[subSpaceNum]->GetTNSEResidual(modeSolution_i, defectMode);
                // correction due to L^2_O Pressure space
                if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
                    IntoL20Vector2D(defectMode + 2 * N_U, N_P, pressure_space_code);

                residual = Ddot(2 * N_U + N_P, defectMode, defectMode);
                impuls_residual = Ddot(2 * N_U + N_P, defectMode, defectMode);
                OutPut("Nonlinear iteration step   0");
                OutPut(setw(14) << impuls_residual);
                OutPut(setw(14) << residual - impuls_residual);
                OutPut(setw(14) << sqrt(residual) << endl);

                //======================================================================
                // Solve the system
                // Nonlinear iteration of fixed point type
                //======================================================================
                for (j = 1; j <= Max_It; j++)
                {
                    // Solve the NSE system
                    SystemMatrixModeAll[subSpaceNum]->Solve(modeSolution_i);

                    if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
                        IntoL20FEFunction(modeSolution_i + 2 * N_U, N_P, Pressure_FeSpace, velocity_space_code, pressure_space_code);

                    // no nonlinear iteration for Stokes problem
                    if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == STOKES)
                        break;

                    // restore the mass matrix for the next nonlinear iteration
                    SystemMatrixModeAll[subSpaceNum]->RestoreMassMat();

                    // assemble the system matrix with given aux, sol and rhs
                    SystemMatrixModeAll[subSpaceNum]->AssembleANonLinear(modeSolution_i, modeSolution_rhs);

                    // assemble system mat, S = M + dt\theta_1*A
                    SystemMatrixModeAll[subSpaceNum]->AssembleSystMatNonLinear();

                    // get the residual
                    // calculate the residual
                    memset(defectMode, 0, (2 * N_U + N_P) * SizeOfDouble);

                    SystemMatrixModeAll[subSpaceNum]->GetTNSEResidual(modeSolution_i, defectMode);

                    // correction due to L^2_O Pressure space
                    if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
                        IntoL20Vector2D(defectMode + 2 * N_U, N_P, pressure_space_code);

                    residual = Ddot(2 * N_U + N_P, defectMode, defectMode);
                    impuls_residual = Ddot(2 * N_U + N_P, defectMode, defectMode);
                    OutPut("Nonlinear iteration step   0");
                    OutPut(setw(14) << impuls_residual);
                    OutPut(setw(14) << residual - impuls_residual);
                    OutPut(setw(14) << sqrt(residual) << endl);

                    if (sqrt(residual) <= limit)
                        break;

                } // for(j=1;j<=Max_It;j++)

                SystemMatrixModeAll[subSpaceNum]->RestoreMassMat();

            } // Subspace Loop Ends

            for (int j = 0; j < subDim; j++)
            {
                for (int i = 0; i < N_U; i++)
                {
                    solMode1[j * N_U + i] = solModeAll[(2 * j * N_U) + i];
                }
            }
            reorthonormalizeC(solMode1, N_U, subDim);
            for (int j = 0; j < subDim; j++)
            {
                for (int i = 0; i < N_U; i++)
                {
                    solModeAll[(2 * j * N_U) + i] = solMode1[j * N_U + i];
                }
            }
            memcpy(solMode, solModeAll, N_Total_ModeDOF * SizeOfDouble);
            // for (int i = 0; i < N_U; i++)
            // {
            // 	for (int j = 0; j < subDim; j++)
            // 	{
            // 		u1_Mode[j * N_U + i] = solModeAll[(2 * j * N_U) + i]; //??
            // 		u2_Mode[j * N_U + i] = 0;							  //??
            // 	}
            // }

            // calcIPMatx(IPMatxMode, u1_Mode, N_U, subDim, 'C');
            // fileoutIPMode = generateFileName(IPModeBaseName, m, N_Realisations);
            // printToTxt(fileoutIPMode, IPMatxMode, subDim, subDim, 'C');

        } // for(l=0;l<N_SubSteps;

        //======================================================================
        // produce outout
        //======================================================================
        if (m == 1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
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

        // start here 12-05-21
        // xxxxxxxxxxxxxxxxxxxxxxxxxxxx
        // fileoutCoeff = generateFileName(coeffBaseName, m, N_Realisations);
        // printToTxt(fileoutCoeff, CoeffVector, N_Realisations, subDim, 'C');

        //======================================================================
        // produce outout
        //======================================================================
        if (m == 1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
        {
            for (int sd = 0; sd < subDim; sd++)
            {

                std::string filenameMode = "Mode_NRealisations_" + std::to_string(N_Realisations) + "ModeN0_" + std::to_string(sd);
                VtkBaseNameMode = const_cast<char *>(filenameMode.c_str());
                if (TDatabase::ParamDB->WRITE_VTK)
                {
                    os.seekp(std::ios::beg);
                    if (modeimg[sd] < 10)
                        os << "VTK/" << VtkBaseNameMode << ".0000" << modeimg[sd] << ".vtk" << ends;
                    else if (modeimg[sd] < 100)
                        os << "VTK/" << VtkBaseNameMode << ".000" << modeimg[sd] << ".vtk" << ends;
                    else if (modeimg[sd] < 1000)
                        os << "VTK/" << VtkBaseNameMode << ".00" << modeimg[sd] << ".vtk" << ends;
                    else if (modeimg[sd] < 10000)
                        os << "VTK/" << VtkBaseNameMode << ".0" << modeimg[sd] << ".vtk" << ends;
                    else
                        os << "VTK/" << VtkBaseNameMode << "." << modeimg[sd] << ".vtk" << ends;
                    OutputModeAll[sd]->WriteVtk(os.str().c_str());
                    modeimg[sd]++;
                }
            }
        }

    } // while(TDatabase::TimeDB->CURRENTTIME< e)

    //======================================================================
    // produce final outout
    //======================================================================
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

    // for (int sd = 0; sd < subDim; sd++)
    // {

    // 	std::string filenameMode = "Mode_NRealisations_" + std::to_string(N_Realisations) + "_ModeN0_" + std::to_string(sd);
    // 	VtkBaseNameMode = const_cast<char *>(filenameMode.c_str());
    // 	if (TDatabase::ParamDB->WRITE_VTK)
    // 	{
    // 		os.seekp(std::ios::beg);
    // 		if (modeimg[sd] < 10)
    // 			os << "VTK/" << VtkBaseNameMode << ".0000" << modeimg[sd] << ".vtk" << ends;
    // 		else if (modeimg[sd] < 100)
    // 			os << "VTK/" << VtkBaseNameMode << ".000" << modeimg[sd] << ".vtk" << ends;
    // 		else if (modeimg[sd] < 1000)
    // 			os << "VTK/" << VtkBaseNameMode << ".00" << modeimg[sd] << ".vtk" << ends;
    // 		else if (modeimg[sd] < 10000)
    // 			os << "VTK/" << VtkBaseNameMode << ".0" << modeimg[sd] << ".vtk" << ends;
    // 		else
    // 			os << "VTK/" << VtkBaseNameMode << "." << modeimg[sd] << ".vtk" << ends;
    // 		OutputModeAll[sd]->WriteVtk(os.str().c_str());
    // 		modeimg[sd]++;
    // 	}
    // }

    TDatabase::TimeDB->CURRENTTIME = 0;

    CloseFiles();

    return 0;
}