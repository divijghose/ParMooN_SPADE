// =======================================================================
//
// Purpose:     main program for solving a time-dependent NSE equation in ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 03.09.2014

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
// include current example
// =======================================================================
#include "../Examples/TNSE_2D/DrivenCavity.h" //   in unit square
// #include "../Examples/TNSE_2D/Bsp3.h" // smooth sol in unit square
// #include "../Examples_All/TNSE_2D/Benchmark2.h"
//#include "../Examples/TNSE_2D/SinCos.h" // smooth sol in unit square
//#include "../Examples/TNSE_2D/SinCos2.h" // smooth sol in unit square
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
    int N_Realisations  = TDatabase::ParamDB->REALIZATIONS;
    double LengthScale  = TDatabase::ParamDB->LENGTHSCALE;
    double EigenPercent = TDatabase::ParamDB->EIGENPERCENT;

	double *org_x_coord = new double[N_U];
	double *org_y_coord = new double[N_U];
	double *x_coord = new double[N_U];
	double *y_coord = new double[N_U];
	int *mappingArray = new int[N_U];

	i = 0;
	int N = (2 * pow(2, TDatabase::ParamDB->UNIFORM_STEPS)) + 1;
	for (int i = 0; i < N_U; i++)
	{
		int local_i = i / N;
		int local_j = i % N;

		x_coord[i] = double(1.0 / (N - 1)) * local_i;
		y_coord[i] = double(1.0 / (N - 1)) * local_j;
	}
	cout << " End File Read" << endl;
	Velocity_FeSpace->GetDOFPosition(org_x_coord, org_y_coord);

	for (int i = 0; i < N_U; i++) // Generated Values
	{
		// get the generated Value
		double xx = x_coord[i];
		double yy = y_coord[i];
		bool foundFlag = false;

		for (int j = 0; j < N_U; j++) // Actual parmooN Co-ordinates
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

	double *x = new double[N_U];
	double *y = new double[N_U];

	for (int i = 0; i < N_U; i++)
	{
		int local_i = i / N;
		int local_j = i % N;

		x[i] = double(1.0 / (N - 1)) * local_j;
		y[i] = double(1.0 / (N - 1)) * local_i;
	}

	double *C = new double[N_U * N_U];	//MATRIX
	double *C1 = new double[N_U * N_U]; //MATRIX  - Corelation Matrix
	double norm = 0;
	for (int i = 0; i < N_U; i++)
	{
		double actual_x = x[i];
		double actual_y = y[i];

		for (int j = 0; j < N_U; j++)
		{
			double local_x = x[j];
			double local_y = y[j];

			double r = sqrt(pow((actual_x - local_x), 2) + pow((actual_y - local_y), 2));

			// CO -Relation
			C[j * N_U + i] = exp((-1.0 * r) / (LengthScale));
			C1[j * N_U + i] = exp((-1.0 * r) / (LengthScale));

		if(TDatabase::ParamDB->stddev_switch == 0)
            {
                double sig_r1 = exp (-1.0/(1.0 - pow(( 2*actual_x - 1),4) ) )  * exp ( -1.0/ ( 1 - pow(( 2*actual_y - 1),4) ) ) ;
                double sig_r2 = exp (-1.0/(1.0 - pow(( 2*local_x - 1),4) ) )  * exp ( -1.0/ ( 1 - pow(( 2*local_y - 1),4) ) ) ; 
            
                // Co Variance
                C[j*N_U + i] *= sig_r1 * sig_r2 * 5.0;
            }

            else if(TDatabase::ParamDB->stddev_switch == 1)
            {
                double E = TDatabase::ParamDB->stddev_denom;
                double disp = TDatabase::ParamDB->stddev_disp;
                double power = TDatabase::ParamDB->stddev_power;
                double sig_r1 = exp ( - pow( ( 2*actual_x - 1 - disp),power)  / (E) )  / (2*3.14159265359 * sqrt(E))  * exp ( - pow(( 2*actual_x - 1-disp),power)  / (E) )  / (2*3.14159265359 * sqrt(E)) ;
                double sig_r2 = exp ( - pow(( 2*local_x - 1 -disp),power)  / (E) )  / (2*3.14159265359 * sqrt(E))  * exp ( - pow(( 2*local_y - 1-disp),power)  / (E) ) / (2*3.14159265359 * sqrt(E)); 
                // Co Variance
                C[j*N_U + i] *= sig_r1 * sig_r2 ;
            }

            else{
                cout << "Error " <<endl;
                exit(0);
            }

            norm += C[j*N + i]*C[j*N + i];
        }

    }

	std::ofstream fileo;
	fileo.open("Corelation.txt");

	for (int i = 0; i < N_U; i++)
	{
		for (int j = 0; j < N_U; j++)
		{
			fileo << C1[i * N_U + j];
			if (j != N_U - 1)
				fileo << ",";
		}
		fileo << endl;
	}

	fileo.close();

	std::ofstream fileo_r;
	fileo.open("Covarriance.txt");

	for (int i = 0; i < N_U; i++)
	{
		for (int j = 0; j < N_U; j++)
		{
			fileo_r << C[i * N_U + j];
			if (j != N_U - 1)
				fileo_r << ",";
		}
		fileo_r << endl;
	}

	fileo_r.close();

	/////////////////////////////SVD//////////////////////////////////////////////
	// Declare SVD parameters
	MKL_INT m1 = N_U, n = N_U, lda = N_U, ldu = N_U, ldvt = N_U, info;
	double superb[std::min(N_U, N_U) - 1];

	double *S = new double[N_U];
	double *U = new double[N_U * N_U];
	double *Vt = new double[N_U * N_U];
	cout << " REALISATIONS COMPUTED " << endl;
	info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', m1, n, C, lda,
						  S, U, ldu, Vt, ldvt, superb);

	cout << endl
		 << endl;

	if (info > 0)
	{
		printf("The algorithm computing SVD failed to converge.\n");
		exit(1);
	}
	cout << " REALISATIONS COMPUTED " << endl;
	int energyVal = 0;

	double sumSingularVal = 0;
	for (int i = 0; i < N_U; i++)
		sumSingularVal += S[i];
	double val = 0;
	for (energyVal = 0; energyVal < N_U; energyVal++)
	{
		val += S[energyVal];
		if (val / sumSingularVal > 0.8)
			break;
	}

	cout << " MODES : " << energyVal + 1 << endl;

	int modDim = energyVal + 1;

	double *Ut = new double[N_U * modDim]();
	double *Z = new double[N_Realisations * modDim]();

	double *SolutionVector = new double[N_U * N_Realisations]();
	// -------------- Generate Random Number Based on Normal Distribution -------------------------//
	int k = 0;
	int skip = N_U - modDim;
	int count = 0;

	for (int i = 0; i < N_U * N_U; i++)
	{
		// cout << "i val " << i <<endl;
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

	///////////
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
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_U, N_Realisations, modDim, 1.0, Ut, modDim, Z, N_Realisations, 0.0, SolutionVector, N_Realisations);
	cout << " MULT DONE " << endl;
	// printMatrix(SolutionVector, N_DOF,N_Realisations);

	// mkl_dimatcopy('R','T', N_DOF,N_Realisations,1.0,SolutionVector,N_DOF,N_Realisations);
	cout << " COPY DONE " << endl;

	cout << " REALISATIONS COMPUTED " << endl;

	//////////////////////////////////End of Realization/////////////////////////////////////////

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

	//define the aux
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
			sol[mappingArray[i]] = SolutionVector[RealNo + N_Realisations * i];
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

				//copy sol, rhs to olssol, oldrhs
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

				//scale B matices and assemble NSE-rhs based on the \theta time stepping scheme
				SystemMatrix->AssembleSystMat(tau / oldtau, oldrhs, rhs, sol);
				oldtau = tau;

				// calculate the residual
				defect = new double[N_TotalDOF];
				memset(defect, 0, N_TotalDOF * SizeOfDouble);

				SystemMatrix->GetTNSEResidual(sol, defect);

				//correction due to L^2_O Pressure space
				if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
					IntoL20Vector2D(defect + 2 * N_U, N_P, pressure_space_code);

				residual = Ddot(N_TotalDOF, defect, defect);
				impuls_residual = Ddot(2 * N_U, defect, defect);
				OutPut("Nonlinear iteration step   0");
				OutPut(setw(14) << impuls_residual);
				OutPut(setw(14) << residual - impuls_residual);
				OutPut(setw(14) << sqrt(residual) << endl);

				//======================================================================
				//Solve the system
				//Nonlinear iteration of fixed point type
				//======================================================================
				for (j = 1; j <= Max_It; j++)
				{
					// Solve the NSE system
					SystemMatrix->Solve(sol);

					if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
						IntoL20FEFunction(sol + 2 * N_U, N_P, Pressure_FeSpace, velocity_space_code, pressure_space_code);

					//no nonlinear iteration for Stokes problem
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

					//correction due to L^2_O Pressure space
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

				//copy sol, rhs to olssol, oldrhs
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

				//scale B matices and assemble NSE-rhs based on the \theta time stepping scheme
				SystemMatrix->AssembleSystMat(tau / oldtau, oldrhs, rhs, sol);
				oldtau = tau;

				// calculate the residual
				defect = new double[N_TotalDOF];
				memset(defect, 0, N_TotalDOF * SizeOfDouble);

				SystemMatrix->GetTNSEResidual(sol, defect);

				//correction due to L^2_O Pressure space
				if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
					IntoL20Vector2D(defect + 2 * N_U, N_P, pressure_space_code);

				residual = Ddot(N_TotalDOF, defect, defect);
				impuls_residual = Ddot(2 * N_U, defect, defect);
				OutPut("Nonlinear iteration step   0");
				OutPut(setw(14) << impuls_residual);
				OutPut(setw(14) << residual - impuls_residual);
				OutPut(setw(14) << sqrt(residual) << endl);

				//======================================================================
				//Solve the system
				//Nonlinear iteration of fixed point type
				//======================================================================
				for (j = 1; j <= Max_It; j++)
				{
					// Solve the NSE system
					SystemMatrix->Solve(sol);

					if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
						IntoL20FEFunction(sol + 2 * N_U, N_P, Pressure_FeSpace, velocity_space_code, pressure_space_code);

					//no nonlinear iteration for Stokes problem
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

					//correction due to L^2_O Pressure space
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
	pyRead << N_U << endl << N_P << endl << N_Realisations << endl << m;

	CloseFiles();
	return 0;
} // end main
