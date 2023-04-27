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

// #include "../Examples/TNSE_2D/Bsp3.h" // smooth sol in unit square
// #include "../Examples_All/TNSE_2D/Benchmark2.h"
// #include "../Examples/TNSE_2D/SinCos.h" // smooth sol in unit square
// #include "../Examples/TNSE_2D/SinCos2.h" // smooth sol in unit square
// =======================================================================
// main program
// =======================================================================
int main(int argc, char *argv[])
{
    // ======================================================================
    //  declaration of variables
    // ======================================================================
    int i, j, l, m, N_Cells, ORDER, N_U, N_P, N_L, N_TotalDOF, img = 1, pressure_space_code;
    int Max_It, NSEType, velocity_space_code, N_SubSteps, Disctype;

    double *sol, *rhs, *oldrhs, *defect, t1, t2, residual, impuls_residual;
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
    char *PsBaseName, *VtkBaseName, *GEO;
    char UString[] = "u";
    char PString[] = "p";
    char NameString[] = "VMS";

    std::ostringstream os;
    os << " ";

    mkdir(vtkdir, 0777);

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

    //======================================================================
    // construct all finite element functions
    //======================================================================
    sol = new double[N_TotalDOF];
    rhs = new double[N_TotalDOF];
    oldrhs = new double[N_TotalDOF];

    memset(sol, 0, N_TotalDOF * SizeOfDouble);
    memset(rhs, 0, N_TotalDOF * SizeOfDouble);

    Velocity = new TFEVectFunct2D(Velocity_FeSpace, UString, UString, sol, N_U, 2);
    u1 = Velocity->GetComponent(0);
    u2 = Velocity->GetComponent(1);
    Pressure = new TFEFunction2D(Pressure_FeSpace, PString, PString, sol + 2 * N_U, N_P);

    //  interpolate the initial solution
    u1->Interpolate(InitialU1);
    u2->Interpolate(InitialU2);
    Pressure->Interpolate(InitialP);

    ///////////////////////////////////////////////////////////////////////////////////////////////
    ////////// -------- REALISATION DATA GENERATION ----------------------------------------- //////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    int N_Realisations = TDatabase::ParamDB->REALIZATIONS;
    double *RealizationVector = new double[N_U * N_Realisations]();

    GenerateRealizations(Velocity_FeSpace, RealizationVector);
    
    //////////////////////////////////End of Realization/////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////////////////
    /////////// DO - Initialization /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////

    double *MeanVector = new double[N_U * 1]();
    int subDim = calculateStochSubspaceDim(Velocity_FeSpace, RealizationVector);

    double *CoeffVector = new double[N_Realisations * subDim]();
    double *ModeVector = new double[N_U * subDim]();
    InitializeDO(Velocity_FeSpace, RealizationVector, MeanVector, ModeVector, CoeffVector);
    cout << "Works till Initialize" << endl;
    
    ////////////////////////////////////////////DO - Initialization Ends//////////////////////////////////////
    ///////================================================================================//////////////////

    CalcCovarianceMatx(CoeffVector); // Added to Database.h

    //=========================================================================
    // Assign dimension values to Database
    //=========================================================================
    TDatabase::ParamDB->N_Subspace_Dim = subDim; // Added to Database.h
    TDatabase::ParamDB->REALIZATIONS = N_Realisations;

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
    double *TotalModeVector = new double[N_U * N_Realisations * 2]();

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
    // ******************** Solve Mean Equation ********************
    //======================================================================

    //======================================================================
    // construct all finite element functions
    //======================================================================
    double *solMean, *rhsMean, *oldrhsMean; //,*defect, t1, t2, residual, impuls_residual;
    solMean = new double[N_TotalDOF];
    rhsMean = new double[N_TotalDOF];
    oldrhsMean = new double[N_TotalDOF];

    memset(solMean, 0, N_TotalDOF * SizeOfDouble);
    memset(rhsMean, 0, N_TotalDOF * SizeOfDouble);

    // TFESpace2D *Velocity_FeSpace, *Pressure_FeSpace, *fesp[2];

    // Initialize New VectFunctions and FEFunctions for Mean Vector Solution
    TFEVectFunct2D *VelocityMean;
    TFEFunction2D *u1Mean, *u2Mean, *PressureMean;
    // TOutput2D *Output;
    TSystemTNSE2D *SystemMatrixMean;
    // TAuxParam2D *aux, *NSEaux_error;
    // MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
    VelocityMean = new TFEVectFunct2D(Velocity_FeSpace, UString, UString, solMean, N_U, 2);
    u1Mean = VelocityMean->GetComponent(0);
    u2Mean = VelocityMean->GetComponent(1);
    PressureMean = new TFEFunction2D(Pressure_FeSpace, PString, PString, sol + 2 * N_U, N_P);

    // Interpolate the initial solution
    u1Mean->Interpolate(InitialU1);
    u2Mean->Interpolate(InitialU2);
    PressureMean->Interpolate(InitialP);

    // ======================================================================
    // SystemMatrix construction and solution
    // ======================================================================

    SystemMatrixMean = new TSystemTNSE2D(Velocity_FeSpace, Pressure_FeSpace, VelocityMean, PressureMean, solMean, rhsMean, Disctype, NSEType, DIRECT
#ifdef __PRIVATE__
                                         ,
                                         Projection_space, NULL, NULL
#endif
    );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //======================================================================
    // SystemMatrix construction and solution
    //======================================================================
    // Disc type: GALERKIN
    // Solver: AMG_SOLVE (or) GMG  (or) DIRECT
    SystemMatrix = new TSystemTNSE2D(Velocity_FeSpace, Pressure_FeSpace, Velocity, Pressure, sol, rhs, Disctype, NSEType, DIRECT
#ifdef __PRIVATE__
                                     ,
                                     Projection_space, NULL, NULL
#endif
    );

    // define the aux
    fesp[0] = Velocity_FeSpace;

    fefct[0] = u1;
    fefct[1] = u2;

    switch (Disctype) /// Shouldn't there be a switch case for Disctype = DO_Disctype
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
    OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
    //      TDatabase::TimeDB->TIMESTEPLENGTH = hmin;
    //       cout<<TDatabase::TimeDB->TIMESTEPLENGTH<<"\n";

    //======================================================================

    // initilize the system matrix with the functions defined in Example file
    // last argument is aux that is used to pass additional fe functions (eg. mesh velocity)
    SystemMatrix->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue, aux, NSEaux_error);
    /////////////////////////////////////////Monte Carlo//////////////////////////////////////

    // Setup array for random number
    srand(time(NULL));
    int N_samples = 100;
    int *indexArray = new int[N_samples];
    for (int i = 0; i < N_samples; i++)
        indexArray[i] = rand() % N_U;

    ///////////////////////////Looping for Realization////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    std::ofstream fileoutSolution;
    fileoutSolution.open("FinalSolution_AllRealisations.txt");
    std::ofstream fileout;

    for (int RealNo = 0; RealNo < N_Realisations; RealNo++)
    {
        cout << " Real no " << RealNo << endl;
        ////////////////////////Divergence Free Adjustment - Run for one time step//////////////////////
        // assemble M, A matrices and rhs

        for (int i = 0; i < N_U; i++)
            sol[i] = RealizationVector[RealNo + N_Realisations * i];
        SystemMatrix->Assemble(sol, rhs);

        VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
        std::string str = std::to_string(RealNo);
        std::string filename = "Realization_Nr_" + std::to_string(RealNo);
        VtkBaseName = const_cast<char *>(filename.c_str());

        std::string name = "Realization _Number_" + std::to_string(int(RealNo));
        mkdir(name.c_str(), 0777);

        //======================================================================
        // time disc loop
        //======================================================================
        // parameters for time stepping scheme
        m = 0;
        N_SubSteps = 2;
        oldtau = 1.;
        end_time = 0.01;
        limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
        Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
        memset(AllErrors, 0, 7. * SizeOfDouble);

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

                tau = 0.005;
                TDatabase::TimeDB->CURRENTTIME += tau;

                OutPut(endl
                       << "CURRENT TIME: ");
                OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

                // copy sol, rhs to olssol, oldrhs
                memcpy(oldrhs, rhs, N_TotalDOF * SizeOfDouble);

                // assemble only rhs, nonlinear matrix for NSE will be assemble in fixed point iteration
                // not needed if rhs is not time-dependent
                if (m != 1)
                {
                    SystemMatrix->AssembleRhs(sol, rhs);
                }
                else
                {
                    SystemMatrix->Assemble(sol, rhs);
                }

                // scale B matices and assemble NSE-rhs based on the \theta time stepping scheme
                SystemMatrix->AssembleSystMat(tau / oldtau, oldrhs, rhs, sol);
                oldtau = tau;

                // calculate the residual
                defect = new double[N_TotalDOF];
                memset(defect, 0, N_TotalDOF * SizeOfDouble);

                SystemMatrix->GetTNSEResidual(sol, defect);

                // correction due to L^2_O Pressure space
                if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
                    IntoL20Vector2D(defect + 2 * N_U, N_P, pressure_space_code);

                residual = Ddot(N_TotalDOF, defect, defect);
                impuls_residual = Ddot(2 * N_U, defect, defect);
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
                    SystemMatrix->Solve(sol);

                    if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
                        IntoL20FEFunction(sol + 2 * N_U, N_P, Pressure_FeSpace, velocity_space_code, pressure_space_code);

                    // no nonlinear iteration for Stokes problem
                    if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == STOKES)
                        break;

                    // restore the mass matrix for the next nonlinear iteration
                    SystemMatrix->RestoreMassMat();

                    // assemble the system matrix with given aux, sol and rhs
                    SystemMatrix->AssembleANonLinear(sol, rhs);

                    // assemble system mat, S = M + dt\theta_1*A
                    SystemMatrix->AssembleSystMatNonLinear();

                    // get the residual
                    memset(defect, 0, N_TotalDOF * SizeOfDouble);
                    SystemMatrix->GetTNSEResidual(sol, defect);

                    // correction due to L^2_O Pressure space
                    if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
                        IntoL20Vector2D(defect + 2 * N_U, N_P, pressure_space_code);

                    residual = Ddot(N_TotalDOF, defect, defect);
                    impuls_residual = Ddot(2 * N_U, defect, defect);
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
                SystemMatrix->RestoreMassMat();

            } // for(l=0;l<N_SubSteps;
              //======================================================================
              // measure errors to known solution
              //======================================================================
            if (TDatabase::ParamDB->MEASURE_ERRORS)
            {
                //        u1->Interpolate(ExactU1);
                //        u2->Interpolate(ExactU2);
                //        Pressure->Interpolate(ExactP);

                SystemMatrix->MeasureTNSEErrors(ExactU1, ExactU2, ExactP, AllErrors);

                OutPut("L2(u): " << AllErrors[0] << endl);
                OutPut("H1-semi(u): " << AllErrors[1] << endl);
                OutPut("L2(p): " << AllErrors[2] << endl);
                OutPut("H1-semi(p): " << AllErrors[3] << endl);
                OutPut(AllErrors[4] << " l_infty(L2(u)) " << AllErrors[5] << endl);
                OutPut(TDatabase::TimeDB->CURRENTTIME << " L2(0,t,L2)(u) : " << sqrt(AllErrors[6]) << endl);

            } // if(TDatabase::ParamDB->MEASURE_ERRORS)

            //======================================================================
            // produce outout
            //======================================================================

        } // while(TDatabase::TimeDB->CURRENTTIME< e

        //======================================================================
        // produce final outout
        //======================================================================

        TDatabase::TimeDB->CURRENTTIME = 0.0;
        //////////////////Divergence Adjustment Ends/////////////////////////////////////////////////

        // assemble M, A matrices and rhs
        SystemMatrix->Assemble(sol, rhs);

        //======================================================================
        // produce outout
        //======================================================================

        Output = new TOutput2D(2, 2, 1, 1, Domain);

        Output->AddFEVectFunct(Velocity);
        Output->AddFEFunction(Pressure);

        if (TDatabase::ParamDB->WRITE_VTK)
        {
            os.seekp(std::ios::beg);
            if (img < 10)
                os << name.c_str() << "/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
            else if (img < 100)
                os << "VTK/" << VtkBaseName << ".000" << img << ".vtk" << ends;
            else if (img < 1000)
                os << name.c_str() << "/" << VtkBaseName << ".00" << img << ".vtk" << ends;
            else if (img < 10000)
                os << name.c_str() << "/" << VtkBaseName << ".0" << img << ".vtk" << ends;
            else
                os << name.c_str() << "/" << VtkBaseName << "." << img << ".vtk" << ends;

            Output->WriteVtk(os.str().c_str());

            img++;
        }

        fileout.open(str + ".txt");
        for (int i = 0; i < N_TotalDOF; i++)
            fileout << sol[i] << "\t";
        fileout << endl;
        exit(0);
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
                memcpy(oldrhs, rhs, N_TotalDOF * SizeOfDouble);

                // assemble only rhs, nonlinear matrix for NSE will be assemble in fixed point iteration
                // not needed if rhs is not time-dependent
                if (m != 1)
                {
                    SystemMatrix->AssembleRhs(sol, rhs);
                }
                else
                {
                    SystemMatrix->Assemble(sol, rhs);
                }

                // scale B matices and assemble NSE-rhs based on the \theta time stepping scheme
                SystemMatrix->AssembleSystMat(tau / oldtau, oldrhs, rhs, sol);
                oldtau = tau;

                // calculate the residual
                defect = new double[N_TotalDOF];
                memset(defect, 0, N_TotalDOF * SizeOfDouble);

                SystemMatrix->GetTNSEResidual(sol, defect);

                // correction due to L^2_O Pressure space
                if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
                    IntoL20Vector2D(defect + 2 * N_U, N_P, pressure_space_code);

                residual = Ddot(N_TotalDOF, defect, defect);
                impuls_residual = Ddot(2 * N_U, defect, defect);
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
                    SystemMatrix->Solve(sol);

                    if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
                        IntoL20FEFunction(sol + 2 * N_U, N_P, Pressure_FeSpace, velocity_space_code, pressure_space_code);

                    // no nonlinear iteration for Stokes problem
                    if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == STOKES)
                        break;

                    // restore the mass matrix for the next nonlinear iteration
                    SystemMatrix->RestoreMassMat();

                    // assemble the system matrix with given aux, sol and rhs
                    SystemMatrix->AssembleANonLinear(sol, rhs);

                    // assemble system mat, S = M + dt\theta_1*A
                    SystemMatrix->AssembleSystMatNonLinear();

                    // get the residual
                    memset(defect, 0, N_TotalDOF * SizeOfDouble);
                    SystemMatrix->GetTNSEResidual(sol, defect);

                    // correction due to L^2_O Pressure space
                    if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
                        IntoL20Vector2D(defect + 2 * N_U, N_P, pressure_space_code);

                    residual = Ddot(N_TotalDOF, defect, defect);
                    impuls_residual = Ddot(2 * N_U, defect, defect);
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
                SystemMatrix->RestoreMassMat();

            } // for(l=0;l<N_SubSteps;
              //======================================================================
              // measure errors to known solution
              //======================================================================
            if (TDatabase::ParamDB->MEASURE_ERRORS)
            {
                //        u1->Interpolate(ExactU1);
                //        u2->Interpolate(ExactU2);
                //        Pressure->Interpolate(ExactP);

                SystemMatrix->MeasureTNSEErrors(ExactU1, ExactU2, ExactP, AllErrors);

                OutPut("L2(u): " << AllErrors[0] << endl);
                OutPut("H1-semi(u): " << AllErrors[1] << endl);
                OutPut("L2(p): " << AllErrors[2] << endl);
                OutPut("H1-semi(p): " << AllErrors[3] << endl);
                OutPut(AllErrors[4] << " l_infty(L2(u)) " << AllErrors[5] << endl);
                OutPut(TDatabase::TimeDB->CURRENTTIME << " L2(0,t,L2)(u) : " << sqrt(AllErrors[6]) << endl);

            } // if(TDatabase::ParamDB->MEASURE_ERRORS)

            //======================================================================
            // produce outout
            //======================================================================
            if (m == 1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
                if (TDatabase::ParamDB->WRITE_VTK)
                {
                    os.seekp(std::ios::beg);
                    if (img < 10)
                        os << name.c_str() << "/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
                    else if (img < 100)
                        os << name.c_str() << "/" << VtkBaseName << ".000" << img << ".vtk" << ends;
                    else if (img < 1000)
                        os << name.c_str() << "/" << VtkBaseName << ".00" << img << ".vtk" << ends;
                    else if (img < 10000)
                        os << name.c_str() << "/" << VtkBaseName << ".0" << img << ".vtk" << ends;
                    else
                        os << name.c_str() << "/" << VtkBaseName << "." << img << ".vtk" << ends;

                    Output->WriteVtk(os.str().c_str());
                    img++;
                }

            for (int i = 0; i < N_TotalDOF; i++)
                fileout << sol[i] << "\t";
            fileout << endl;

        } // while(TDatabase::TimeDB->CURRENTTIME< e

        //======================================================================
        // produce final outout
        //======================================================================
        if (TDatabase::ParamDB->WRITE_VTK)
        {
            os.seekp(std::ios::beg);
            if (img < 10)
                os << name.c_str() << "/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
            else if (img < 100)
                os << name.c_str() << "/" << VtkBaseName << ".000" << img << ".vtk" << ends;
            else if (img < 1000)
                os << name.c_str() << "/" << VtkBaseName << ".00" << img << ".vtk" << ends;
            else if (img < 10000)
                os << name.c_str() << "/" << VtkBaseName << ".0" << img << ".vtk" << ends;
            else
                os << name.c_str() << "/" << VtkBaseName << "." << img << ".vtk" << ends;
            Output->WriteVtk(os.str().c_str());
            img++;
        }
        for (int i = 0; i < N_TotalDOF; i++)
            fileout << sol[i] << "\t";
        fileout << endl;
        fileout.close();

        for (int i = 0; i < N_TotalDOF; i++)
            fileoutSolution << sol[i] << ",";
        fileoutSolution << endl;
        TDatabase::TimeDB->CURRENTTIME = 0;
    }
    std::ofstream pyRead;
    pyRead.open("pyReadIn.txt");
    pyRead << N_U << endl
           << N_P << endl
           << N_Realisations << endl
           << m;

    CloseFiles();
    return 0;
} // end main
