/**
 * @file TNSE2D_ParMooN_IVUQ_DO_Test.C
 * @brief Purpose:     Main program for solving the set of dynamically orthogonal field
					   equations for Navier-Stokes' equation in 2D.
					   Features included in this main program:
					   1. Mesh generation and database values assignment
					   2. Construction of finite element spaces and functions
					   3. Generation of Monte Carlo realisations
					   4. Initialization of DO state vector
					   5. System construction and solution

 * @authors Sashikumaar Ganesan
 * @authors Divij Ghose
 * @authors Thivin Anandh
 * @bug No known bugs
 */
// =======================================================================
//
// Purpose:     Main program for solving a time-dependent NSE equation in ParMooN with uncertainty quantification using DO field equations
//
// Author:      Sashikumaar Ganesan, Divij Ghose, Thivin Anandh
//
// History:     Implementation started on 03.09.2014
// 				Implementation of DO equations started on 20.04.2022

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
#include "../Examples/DO_UQ/navierstokes_do_ldc_validation.h" //   in unit square
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
	int i, j, l, m, N_Cells, ORDER, N_U, N_P, N_L, N_TotalDOF, img = 1, pressure_space_code, velocity_space_code;

	int pressure_space_code_mean;
	int pressure_space_code_mode;
	int velocity_space_code_mean;
	int velocity_space_code_mode;

	int Max_It, NSEType, N_SubSteps, Disctype;

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

	char *PsBaseName, *VtkBaseName, *GEO;

	char UString[] = "u";
	char PString[] = "p";
	char NameString[] = "VMS";

	std::ostringstream os;
	os << " ";

	std::string meanBaseName = "Mean/Mean_NRealisations_";
	std::string modeBaseName = "Modes/Mode_NRealisations_";
	std::string coeffBaseName = "Coefficients/Coeff_NRealisations_";
	std::string IPMeanBaseName = "IPMatrices/IPMean/IPMean_NRealisations_";
	std::string IPModeBaseName = "IPMatrices/IPMode/IPMode_NRealisations_";
	std::string fileoutCoeff;
	std::string fileoutMean;
	std::string fileoutMode;
	std::string fileoutMC;
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
	int N_Total_DomainDOF = 2 * N_U + N_P;

	// OutPut("Dof Velocity : " << setw(10) << 2 * N_U << endl);
	// OutPut("Dof Pressure : " << setw(10) << N_P << endl);
	// OutPut("Total Dof all: " << setw(10) << N_TotalDOF << endl);

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

	if (TDatabase::ParamDB->toggleRealznSource == 0)
	{
		double LengthScale = TDatabase::ParamDB->LENGTHSCALE;
		double EigenPercent = TDatabase::ParamDB->EIGENPERCENT;

		double *org_x_coord = new double[N_U]();
		double *org_y_coord = new double[N_U]();
		double *x_coord = new double[N_U]();
		double *y_coord = new double[N_U]();
		int *mappingArray = new int[N_U]();

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
		// int N_U =  N * N;
		double *x = new double[N_U]();
		double *y = new double[N_U]();

		for (int i = 0; i < N_U; i++)
		{
			int local_i = i / N;
			int local_j = i % N;

			x[i] = double(1.0 / (N - 1)) * local_j;
			y[i] = double(1.0 / (N - 1)) * local_i;
		}

		double *C = new double[N_U * N_U]();  // MATRIX
		double *C1 = new double[N_U * N_U](); // MATRIX  - Corelation Matrix

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
				C[i * N_U + j] = exp((-1.0 * r) / (LengthScale)) * (1 + (r / LengthScale) + pow(r / LengthScale, 2));
				C1[i * N_U + j] = exp((-1.0 * r) / (LengthScale)) * (1 + (r / LengthScale) + pow(r / LengthScale, 2));

				if (TDatabase::ParamDB->stddev_switch == 0)
				{
					double sig_r1 = exp(-1.0 / (1.0 - pow((2 * actual_x - 1), 4))) * exp(-1.0 / (1 - pow((2 * actual_y - 1), 4)));
					double sig_r2 = exp(-1.0 / (1.0 - pow((2 * local_x - 1), 4))) * exp(-1.0 / (1 - pow((2 * local_y - 1), 4)));
					// Co Variance
					C[i * N_U + j] *= sig_r1 * sig_r2 * 5.0;
				}

				else if (TDatabase::ParamDB->stddev_switch == 1)
				{
					double E = TDatabase::ParamDB->stddev_denom;
					double disp = TDatabase::ParamDB->stddev_disp;
					double power = TDatabase::ParamDB->stddev_power;
					double sig_r1 = (exp(-1.0 * pow((2 * actual_x - 1 - disp), power) / (E)) / (2 * Pi * sqrt(E))) * (exp(-1.0 * pow((2 * actual_y - 1 - disp), power) / (E)) / (2 * Pi * sqrt(E)));
					double sig_r2 = (exp(-1.0 * pow((2 * local_x - 1 - disp), power) / (E)) / (2 * Pi * sqrt(E))) * (exp(-1.0 * pow((2 * local_y - 1 - disp), power) / (E)) / (2 * Pi * sqrt(E)));
					// Co Variance
					C[i * N_U + j] *= 2 * sig_r1 * sig_r2;
				}

				else if (TDatabase::ParamDB->stddev_switch == 2)
				{
					double amplitude = TDatabase::ParamDB->stddev_power;
					double sig_r1 = (amplitude)*sin(-1.0 * Pi * (2 * actual_x - 2)) * sin(-1.0 * Pi * (2 * actual_y - 2));
					double sig_r2 = (amplitude)*sin(-1.0 * Pi * (2 * local_x - 2)) * sin(-1.0 * Pi * (2 * local_y - 2));
					C[i * N_U + j] *= sig_r1 * sig_r2;
				}

				else
				{
					cout << "Error - No standard deviation function is defined for stddev_switch: " << TDatabase::ParamDB->stddev_switch << endl;
					exit(0);
				}

				// norm += C[j * N + i] * C[j * N + i];
			}
		}

		std::string fileOutCorrelation = "Init/Correlation_" + std::to_string(N_U) + "_NR_" + std::to_string(N_Realisations) + ".txt";
		printToTxt(fileOutCorrelation, C, N_U, N_U, 'R');

		////////////////////////////////////////////////////// SVD ////////////////////////////////////////////
		// Declare SVD parameters
		MKL_INT m1 = N_U, n = N_U, lda = N_U, ldu = N_U, ldvt = N_U, info;
		double superb[std::min(N_U, N_U) - 1];

		double *S = new double[N_U]();
		double *U = new double[N_U * N_U]();
		double *Vt = new double[N_U * N_U]();
		info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', m1, n, C, lda,
							  S, U, ldu, Vt, ldvt, superb);

		if (info > 0)
		{
			printf("The algorithm computing SVD of Covariance Matrix failed to converge.\n");
			exit(1);
		}

		std::string fileOutCorrS = "Init/Corr_S_" + std::to_string(N_U) + "_NR_" + std::to_string(N_Realisations) + ".txt";
		printToTxt(fileOutCorrS, S, N_U, 1, 'C');

		std::string fileOutCorrU = "Init/Corr_U_" + std::to_string(N_U) + "_NR_" + std::to_string(N_Realisations) + ".txt";
		printToTxt(fileOutCorrU, U, N_U, N_U, 'R');

		std::string fileOutCorrVt = "Init/Corr_Vt_" + std::to_string(N_U) + "_NR_" + std::to_string(N_Realisations) + ".txt";
		printToTxt(fileOutCorrVt, Vt, N_U, N_U, 'R');

		int energyVal = 0;
		int temp = 0;

		double sumSingularVal = 0;
		for (int i = 0; i < N_U; i++)
			sumSingularVal += S[i];

		double val = 0;
		for (energyVal = 0; energyVal < N_U; energyVal++)
		{
			val += S[energyVal];
			temp++;
			if (val / sumSingularVal > 0.99)
				break;
		}

		cout << " MODES : " << temp + 1 << endl;

		int modDim = temp + 1;

		double *Ut = new double[N_U * modDim]();
		double *Z = new double[N_Realisations * modDim]();

		double *RealizationVectorTemp = new double[N_U * N_Realisations]();

		// -------------- Generate Random Number Based on Normal Distribution -------------------------//
		int k = 0;
		int skip = N_U - modDim;
		int count = 0;
		for (int i = 0; i < N_U * N_U; i++)
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

			double *norm1 = new double[N_Realisations]();

			for (int n = 0; n < N_Realisations; ++n)
			{
				Z[k * N_Realisations + n] = S[k] * d(gen);
			}
		}

		cout << " N_Realisations : " << N_Realisations << endl;
		cout << " MULT START " << endl;
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_U, N_Realisations, modDim, 1.0, Ut, modDim, Z, N_Realisations, 0.0, RealizationVector, N_Realisations);
		cout << " MULT DONE " << endl;

		for (int i = 0; i < N_U; i++)
		{
			for (int j = 0; j < N_Realisations; j++)
			{
				RealizationVectorTemp[mappingArray[i] * N_Realisations + j] = RealizationVector[j + N_Realisations * i];
			}
		}

		memcpy(RealizationVector, RealizationVectorTemp, N_U * N_Realisations * SizeOfDouble);

		cout << N_Realisations << " REALISATIONS COMPUTED " << endl;

		delete[] Ut;
		delete[] Z;
		delete[] RealizationVectorTemp;
		delete[] S;
		delete[] U;
		delete[] Vt;
		delete[] org_x_coord;
		delete[] org_y_coord;
		delete[] x_coord;
		delete[] y_coord;
		delete[] C;
		delete[] C1;

		if (TDatabase::ParamDB->writeRealznToText == 1)
			writeRealizationToText(RealizationVector, N_Realisations, N_U);
	}
	else if (TDatabase::ParamDB->toggleRealznSource == 1)
		readRealizationFromText(RealizationVector, N_Realisations, N_U);
	else
	{
		cout << "Please select correct value of Realization_Source" << endl
			 << TDatabase::ParamDB->toggleRealznSource << " is not an acceptable value" << endl;
		exit(0);
	}

	/////////////////////////////////////// -------- END OF REALISATION DATA SETS ------------ ////////////////////////////////////////////////////////////////
	if (TDatabase::ParamDB->toggleDivFreeAdj == 1)
	{

		///////////////////----Divergence Free Adjustment - New Routine ------------/////////
		//======================================================================
		// construct all finite element functions
		//======================================================================
		sol = new double[N_Total_DomainDOF]();
		rhs = new double[N_Total_DomainDOF]();
		oldrhs = new double[N_Total_DomainDOF]();
		Velocity = new TFEVectFunct2D(Velocity_FeSpace, UString, UString, sol, N_U, 2);
		u1 = Velocity->GetComponent(0);
		u2 = Velocity->GetComponent(1);
		Pressure = new TFEFunction2D(Pressure_FeSpace, PString, PString, sol + 2 * N_U, N_P);
		//  interpolate the initial solution
		u1->Interpolate(InitialU1);
		u2->Interpolate(InitialU2);
		Pressure->Interpolate(InitialP);

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

		for (int RealNo = 0; RealNo < N_Realisations; RealNo++)
		{ /// Realization Loop Begin
			cout << " Divergence-free adjustment Real no " << RealNo << endl;
			////////////////////////Divergence Free Adjustment - Run for one time step//////////////////////
			// assemble M, A matrices and rhs

			for (int i = 0; i < N_U; i++)
				sol[i] = RealizationVector[RealNo + N_Realisations * i];
			SystemMatrix->Assemble(sol, rhs);

			//======================================================================
			// time disc loop
			//======================================================================
			// parameters for time stepping scheme
			m = 0;
			N_SubSteps = 2;
			oldtau = 1.;
			end_time = TDatabase::TimeDB->DF_ENDTIME;
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
						// OutPut("Theta1: " << TDatabase::TimeDB->THETA1 << endl);
						// OutPut("Theta2: " << TDatabase::TimeDB->THETA2 << endl);
						// OutPut("Theta3: " << TDatabase::TimeDB->THETA3 << endl);
						// OutPut("Theta4: " << TDatabase::TimeDB->THETA4 << endl);
					}

					tau = TDatabase::TimeDB->DF_TIMESTEPLENGTH;
					TDatabase::TimeDB->CURRENTTIME += tau;

					// OutPut(endl
					// 	   << "DF ADJUSTMENT CURRENT TIME: ");
					// OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

					// copy sol, rhs to olssol, oldrhs
					memcpy(oldrhs, rhs, N_Total_DomainDOF * SizeOfDouble);

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
					defect = new double[N_Total_DomainDOF];
					memset(defect, 0, N_Total_DomainDOF * SizeOfDouble);

					SystemMatrix->GetTNSEResidual(sol, defect);

					// correction due to L^2_O Pressure space
					if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
						IntoL20Vector2D(defect + 2 * N_U, N_P, pressure_space_code);

					residual = Ddot(N_Total_DomainDOF, defect, defect);
					impuls_residual = Ddot(2 * N_U, defect, defect);
					// OutPut("Nonlinear iteration step   0");
					// OutPut(setw(14) << impuls_residual);
					// OutPut(setw(14) << residual - impuls_residual);
					// OutPut(setw(14) << sqrt(residual) << endl);

					//======================================================================
					// Solve the system
					// Nonlinear iteration of fixed point type
					//======================================================================
					for (j = 1; j <= Max_It; j++)
					{
						// Solve the NSE system
						SystemMatrix->Solve(sol);
						for (int i = 0; i < N_U; i++)
							RealizationVector[RealNo + N_Realisations * i] = sol[i];

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
						memset(defect, 0, N_Total_DomainDOF * SizeOfDouble);
						SystemMatrix->GetTNSEResidual(sol, defect);

						// correction due to L^2_O Pressure space
						if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
							IntoL20Vector2D(defect + 2 * N_U, N_P, pressure_space_code);

						residual = Ddot(N_Total_DomainDOF, defect, defect);
						impuls_residual = Ddot(2 * N_U, defect, defect);
						// OutPut("nonlinear iteration step " << setw(3) << j);
						// OutPut(setw(14) << impuls_residual);
						// OutPut(setw(14) << residual - impuls_residual);
						// OutPut(setw(14) << sqrt(residual) << endl);

						if (sqrt(residual) <= limit)
							break;

					} // for(j=1;j<=Max_It;j++)
					  /*           cout << " test VHM main " << endl;
	exit(0);      */
					// restore the mass matrix for the next time step
					SystemMatrix->RestoreMassMat();

				} // for(l=0;l<N_SubSteps;

			} // while(TDatabase::TimeDB->CURRENTTIME< e

			//======================================================================
			// produce final outout
			//======================================================================

			TDatabase::TimeDB->CURRENTTIME = 0.0;
			//////////////////Divergence Adjustment Ends/////////////////////////////////////////////////
		} /// Realization Loop End
		  ////------------Divergence Free Adjustment - New Routine Ends -----//////////////
	}
	//

	////////////////////////////////////// -------- START OF DO INITIALIZATION ------------ ////////////////////////////////////////////////////////////////

	double *MeanVector = new double[N_U * 1](); // overline{C}_{dof} = \sum_{i=1}^{N_Realisations}(C^{i}_{dof})/N_Realisations
	for (int i = 0; i < N_U; ++i)
	{
		for (int j = 0; j < N_Realisations; ++j)
		{
			MeanVector[i] += (RealizationVector[i * N_Realisations + j] / N_Realisations);
		}
	}
	printToTxt("Init/Mean.txt", MeanVector, N_U, 1, 'C');

	double *PerturbationVector = new double[N_U * N_Realisations](); // \hat{C}^{i}_{dof} = C^{i}_{dof} - \overline{C}_{dof}
	for (int i = 0; i < N_U; ++i)
	{
		for (int j = 0; j < N_Realisations; ++j)
		{
			PerturbationVector[i * N_Realisations + j] = RealizationVector[i * N_Realisations + j] - MeanVector[i];
		}
	}
	printToTxt("Init/PerturbationVector.txt", PerturbationVector, N_U, N_Realisations, 'R');

	const char vtkdir[] = "VTK";
	const char modedir[] = "Modes";
	const char meandir[] = "Mean";
	const char coeffdir[] = "Coefficients";
	const char mcdir[] = "MonteCarlo";
	const char endir[] = "Energy_Data";
	const char initdir[] = "Init";

	mkdir(vtkdir, 0777);
	mkdir(meandir, 0777);
	mkdir(modedir, 0777);
	mkdir(coeffdir, 0777);
	mkdir(mcdir, 0777);
	mkdir(endir, 0777);
	mkdir(initdir, 0777);
	//================================================================================================
	/////////////////////////////DO - Initialization SVD//////////////////////////////////////////////
	//================================================================================================
	// Declare SVD parameters
	int minDim = std::min(N_U, N_Realisations);
	MKL_INT mDO = N_U, nDO = N_Realisations, ldaDO = N_Realisations, lduDO = minDim, ldvtDO = N_Realisations, infoDO;
	double superbDO[minDim - 1];

	double *PerturbationVectorCopy = new double[N_U * N_Realisations]();
	memcpy(PerturbationVectorCopy, PerturbationVector, N_U * N_Realisations * SizeOfDouble);

	double *Sg = new double[minDim];
	double *L = new double[N_U * minDim];
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

	/// Initialize Coefficient Matrix - First subDim columns of Projection Matrix ////////////////////
	double *CoeffVector = new double[N_Realisations * subDim]();

	////////////Initialize Mode Vector - First subDim columns of Left Singular Vector//////////////////

	double *ModeVector = new double[N_U * subDim]();

	for (int i = 0; i < N_U; i++)
	{
		for (int j = 0; j < subDim; j++)
		{
			ModeVector[j * N_U + i] = L[i * minDim + j]; // ModeVector in Col Major form
		}
	}

	TFESpace2D *VelocityModeDOInit_FeSpace = new TFESpace2D(coll, (char *)"Mode_Init", (char *)"FE Space for Mode Solution", BoundCondition, ORDER, NULL);
	TFEVectFunct2D *Velocity_ModeDOInit = new TFEVectFunct2D(VelocityModeDOInit_FeSpace, (char *)"U_Mode_Init", (char *)"Mode Component", ModeVector, N_U, subDim); // check length ??

	double *normzdModeVector = new double[N_U * subDim]();
	double *tempModeVector = new double[N_U * subDim]();

	normalizeStochasticModes(VelocityModeDOInit_FeSpace, Velocity_ModeDOInit, subDim, tempModeVector);
	for (int i = 0; i < N_U; i++)
	{
		for (int j = 0; j < subDim; j++)
		{
			normzdModeVector[i * subDim + j] = tempModeVector[j * N_U + i];
		}
	}
	double *ProjectionVector = new double[N_Realisations * subDim]();

	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N_Realisations, subDim, N_U, 1.0, PerturbationVector, N_Realisations, normzdModeVector, subDim, 0.0, ProjectionVector, subDim);
	for (int i = 0; i < N_Realisations; i++)
	{
		for (int j = 0; j < subDim; j++)
		{
			CoeffVector[j * N_Realisations + i] = ProjectionVector[i * subDim + j]; // CoeffVector in Col Major form
		}
	}

	memcpy(ModeVector, tempModeVector, N_U * subDim * SizeOfDouble);
	printToTxt("Init/Coeff.txt", CoeffVector, N_Realisations, subDim, 'C');
	printToTxt("Init/Mode.txt", ModeVector, N_U, subDim, 'C');

	////////////////////////////////////////////DO - Initialization Ends//////////////////////////////////////
	///////================================================================================//////////////////
	const char ip[] = "IPMatrices";
	mkdir(ip, 0777);
	const char ipmeandir[] = "IPMatrices/IPMean";
	mkdir(ipmeandir, 0777);
	const char ipmodedir[] = "IPMatrices/IPMode";
	mkdir(ipmodedir, 0777);

	double *IPMatxMode = new double[subDim * subDim]();
	double *IPMatxMean = new double[1 * 1]();
	calcIPMatx(IPMatxMode, ModeVector, N_U, subDim, 'C');
	printToTxt("Init/IPMatxMode_Init.txt", IPMatxMode, subDim, subDim, 'R');

	calcIPMatx(IPMatxMean, MeanVector, N_U, 1, 'C');
	printToTxt("Init/IPMatxMean_Init.txt", IPMatxMean, 1, 1, 'R');
	delete[] PerturbationVector;
	delete[] PerturbationVectorCopy;
	delete[] Sg;
	delete[] L;
	delete[] Rt;
	delete[] ProjectionVector;
	//=========================================================================
	// Assign dimension values to Database
	//=========================================================================
	TDatabase::ParamDB->N_Subspace_Dim = subDim; // Added to Database.h
	TDatabase::ParamDB->REALIZATIONS = N_Realisations;

	//=========================================================================
	// Set up FE Spaces for velocity and pressure
	//=========================================================================
	TFESpace2D *VelocityMean_FeSpace, *PressureMean_FeSpace;
	GetVelocityAndPressureSpace(coll, BoundCondition, mortarcoll, VelocityMean_FeSpace,
								PressureMean_FeSpace, &pressure_space_code_mean,
								TDatabase::ParamDB->VELOCITY_SPACE,
								TDatabase::ParamDB->PRESSURE_SPACE);

	// VelocityMean_FeSpace = new TFESpace2D(coll, (char *)"Mean", (char *)"FE Space for Mean Solution", BoundCondition, ORDER, NULL);

	TFESpace2D *VelocityMode_FeSpace, *PressureMode_FeSpace;
	GetVelocityAndPressureSpace(coll, BoundCondition, mortarcoll, VelocityMode_FeSpace,
								PressureMode_FeSpace, &pressure_space_code_mode,
								TDatabase::ParamDB->VELOCITY_SPACE,
								TDatabase::ParamDB->PRESSURE_SPACE);

	// VelocityMode_FeSpace = new TFESpace2D(coll, (char *)"Mode", (char *)"FE Space for Mode Solution", BoundCondition, ORDER, NULL);

	int N_M = VelocityMode_FeSpace->GetN_DegreesOfFreedom();
	int N_Total_MeanDOF = 2 * N_U + N_P;
	int N_Total_ModeDOF = (2 * N_M + N_P) * subDim;
	N_TotalDOF = N_Total_MeanDOF + N_Total_ModeDOF;
	OutPut("Dof Mean Velocity : " << setw(10) << N_Total_MeanDOF << endl);
	OutPut("Dof Mode Velocity : " << setw(10) << N_Total_ModeDOF << endl);
	OutPut("Total DOF Mean+Mode : " << setw(10) << N_TotalDOF << endl);

	//======================================================================
	// construct all finite element functions
	//======================================================================

	// for Mean
	double *solMean, *rhsMean, *old_rhsMean, *defectMean, *old_solMean;
	solMean = new double[N_Total_MeanDOF]();
	rhsMean = new double[N_Total_MeanDOF]();
	old_rhsMean = new double[N_Total_MeanDOF]();
	old_solMean = new double[N_Total_MeanDOF]();

	for (i = 0; i < N_U; i++)
	{

		solMean[i] = MeanVector[i];
		// solMean[N_U + i] = MeanVector[i];
		// solMean[N_U + i] = 0;
	}

	// for (i = 0; i < N_P; i++)
	// {

	// 	solMean[2 * N_U + i] = 0;
	// }

	TFEVectFunct2D *Velocity_Mean;
	TFEFunction2D *Pressure_Mean;

	Velocity_Mean = new TFEVectFunct2D(VelocityMean_FeSpace, (char *)"U_Mean", (char *)"Mean Component", solMean, N_U, 2); // check length

	Pressure_Mean = new TFEFunction2D(PressureMean_FeSpace, (char *)"P_Mean", (char *)"Mean Component", solMean + 2 * N_U, N_P);

	// Interpolate the initial solution
	//  u1Mean->Interpolate(InitialU1);
	//  u2Mean->Interpolate(InitialU2);
	//  Pressure_Mean->Interpolate(InitialP);

	// for mode
	double *solMode, *rhsMode, *old_rhsMode, *defectMode, *old_solMode;
	solMode = new double[N_Total_ModeDOF]();
	old_solMode = new double[2 * N_M + N_P]();
	rhsMode = new double[N_Total_ModeDOF]();
	old_rhsMode = new double[2 * N_M + N_P]();

	double *solMode1 = new double[N_M * subDim]();
	double *solMode2 = new double[N_M * subDim]();

	double *solModeAll = new double[N_Total_ModeDOF]();

	for (int j = 0; j < subDim; j++)
	{
		for (int i = 0; i < N_M; i++)
		{
			solMode[(j * (2 * N_M + N_P)) + i] = ModeVector[j * N_U + i]; // first component velocity
			solModeAll[(j * (2 * N_M + N_P)) + i] = ModeVector[j * N_U + i];
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
	double *solModeVeloCopy = new double[2 * N_M * subDim]();

	for (int i = 0; i < subDim; i++)
	{
		for (int j = 0; j < 2 * N_M; j++)
			solModeVeloCopy[2 * N_M * i + j] = solMode[i * (2 * N_M + N_P) + j];
	}

	double *solModePressCopy = new double[N_P * subDim]();

	for (int i = 0; i < subDim; i++)
	{
		for (int j = 0; j < N_P; j++)
			solModePressCopy[N_P * i + j] = solMode[i * (2 * N_M + N_P) + 2 * N_M + j];
	}

	TFEVectFunct2D *Velocity_Mode;
	TFEFunction2D *u1Mode, *u2Mode, *Pressure_Mode, *fefctMode[2];

	Velocity_Mode = new TFEVectFunct2D(VelocityMode_FeSpace, (char *)"U_Mode", (char *)"Mode Component", solModeVeloCopy, N_M, 2 * subDim); // check length ??

	// u1Mode = Velocity_Mode->GetComponent(0);
	// u2Mode = Velocity_Mode->GetComponent(1);
	Pressure_Mode = new TFEFunction2D(PressureMode_FeSpace, (char *)"P_Mode", (char *)"Mode Component", solModePressCopy, N_P); // check length ??

	TFEVectFunct2D **VelocityModeAll = new TFEVectFunct2D *[subDim];
	for (int subD = 0; subD < subDim; subD++)
		VelocityModeAll[subD] = new TFEVectFunct2D(VelocityMode_FeSpace, (char *)"U_Mode", (char *)"Mode Component", solModeAll + (subD * (2 * N_M + N_P)), N_M, 2);

	TFEFunction2D **PressureModeAll = new TFEFunction2D *[subDim];
	for (int subD = 0; subD < subDim; subD++)
		PressureModeAll[subD] = new TFEFunction2D(PressureMode_FeSpace, (char *)"P_Mode", (char *)"Mode Component", solModeAll + (subDim * (2 * N_M + N_P)) + 2 * N_M, N_P);

	//======================================================================
	// SystemMatrix construction and solution
	//======================================================================

	// -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
	//------------------------------------------ MEAN EQUATION SETUP -----------------------------------------------------//
	TSystemTNSE2D *SystemMatrix_Mean;
	SystemMatrix_Mean = new TSystemTNSE2D(VelocityMean_FeSpace, PressureMean_FeSpace, Velocity_Mean, Pressure_Mean, solMean, rhsMean, Disctype, NSEType, DIRECT
#ifdef __PRIVATE__
										  ,
										  Projection_space, NULL, NULL
#endif
	);

	TFESpace2D *fespMean[2];

	fespMean[0] = VelocityMean_FeSpace;

	TFEFunction2D *u1Mean, *u2Mean, *fefctMean[2];
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
		SystemMatrixModeAll[subD] = new TSystemTNSE2D(VelocityMode_FeSpace, PressureMode_FeSpace, VelocityModeAll[subD], PressureModeAll[subD], solModeAll + (subD * (2 * N_M + N_P)), rhsModeAll + (subD * (2 * N_M + N_P)), Disctype, NSEType, DIRECT
#ifdef __PRIVATE__
													  ,
													  Projection_space, NULL, NULL
#endif
		);

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

	// -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
	//------------------------------------------ MODE EQUATION SETUP END-----------------------------------------------------//
	// -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
	//------------------------------------------ COEFF EQUATION SETUP END-----------------------------------------------------//

	TFEVectFunct2D *FeVector_Coefficient = new TFEVectFunct2D(Velocity_FeSpace, (char *)"Coeff", (char *)"Coefficients", CoeffVector, N_Realisations, subDim);
	// -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
	//------------------------------------------ COEFF EQUATION SETUP END END-----------------------------------------------------//

	// -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
	//------------------------------------------ MEAN EQUATION INITIALIZATION -----------------------------------------------------//
	SystemMatrix_Mean->Init(DO_Mean_Equation_Coefficients, BoundConditionMean, U1BoundValueMean, U2BoundValueMean, auxMean, NSEaux_error_mean);
	SystemMatrix_Mean->Assemble(solMean, rhsMean);

	// -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
	//------------------------------------------ MODE EQUATION INITIALIZATION -----------------------------------------------------//
	for (int subD = 0; subD < subDim; subD++)
	{
		SystemMatrixModeAll[subD]->Init(DO_Mode_Equation_Coefficients, BoundCondition, U1BoundValue, U2BoundValue, auxModeAll[subD], NSEaux_error_modeAll[subD]);
		SystemMatrixModeAll[subD]->Assemble(solModeAll + (subD * (2 * N_M + N_P)), rhsModeAll + (subD * (2 * N_M + N_P))); // seg fault
	}
	TDatabase::ParamDB->COVARIANCE_MATRIX_DO = new double[subDim * subDim]();
	TDatabase::ParamDB->COSKEWNESS_MATRIX_DO = new double[subDim * subDim * subDim]();
	TDatabase::ParamDB->COVARIANCE_INVERSE_DO = new double[subDim * subDim]();

	CalcCovarianceMatx(CoeffVector);
	CalcCoskewnessMatx(CoeffVector);
	InvertCov();

	TDatabase::TimeDB->CURRENTTIME = 0.0;

	// done here
	SystemMatrix_Mean->Assemble(solMean, rhsMean);
	for (int subD = 0; subD < subDim; subD++)
		SystemMatrixModeAll[subD]->Assemble(solModeAll + (subD * (2 * N_M + N_P)), rhsModeAll + (subD * (2 * N_M + N_P)));

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
	// OutputMean->AddFEFunction(Pressure_Mean);

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
			DO_Mean_RHS(VelocityMean_FeSpace, Velocity_Mode, subDim, rhsMean, N_U);

			SystemMatrix_Mean->Assemble(solMean, rhsMean);

			// scale B matices and assemble NSE-rhs based on the \theta time stepping scheme
			SystemMatrix_Mean->AssembleSystMat(tau / oldtau, old_rhsMean, rhsMean, solMean);
			oldtau = tau;

			// calculate the residual

			memset(defectMean, 0, N_Total_MeanDOF * SizeOfDouble);

			SystemMatrix_Mean->GetTNSEResidual(solMean, defectMean);

			// correction due to L^2_O Pressure space
			if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
				IntoL20Vector2D(defectMean + 2 * N_U, N_P, pressure_space_code_mean);

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
					IntoL20FEFunction(solMean + 2 * N_U, N_P, PressureMean_FeSpace, velocity_space_code_mean, pressure_space_code_mean);

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
					IntoL20Vector2D(defectMean + 2 * N_U, N_P, pressure_space_code_mean);

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
				double *modeSolution_i = solMode + subSpaceNum * (2 * N_M + N_P); // This works for column major
				double *modeSolution_rhs = rhsMode + subSpaceNum * (2 * N_M + N_P);
				// copy sol, rhs to olssol, oldrhs
				memcpy(old_rhsMode, modeSolution_rhs, (2 * N_M + N_P) * SizeOfDouble);
				// memcpy(old_solMode, modeSolution_i, N_Total_MeanDOF * SizeOfDouble);

				// Assemble RHS
				DO_Mode_RHS(VelocityMode_FeSpace, Velocity_Mean, Velocity_Mode, Pressure_Mode, subDim, modeSolution_rhs, subSpaceNum);

				SystemMatrixModeAll[subSpaceNum]->Assemble(modeSolution_i, modeSolution_rhs);
				// scale B matices and assemble NSE-rhs based on the \theta time stepping scheme

				SystemMatrixModeAll[subSpaceNum]->AssembleSystMat(tau / oldtau, old_rhsMode, modeSolution_rhs, modeSolution_i);
				oldtau = tau;

				// calculate the residual
				defectMode = new double[2 * N_M + N_P]();

				memset(defectMode, 0, (2 * N_M + N_P) * SizeOfDouble);

				SystemMatrixModeAll[subSpaceNum]->GetTNSEResidual(modeSolution_i, defectMode);
				// correction due to L^2_O Pressure space
				if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
					IntoL20Vector2D(defectMode + 2 * N_M, N_P, pressure_space_code_mode);

				residual = Ddot(2 * N_M + N_P, defectMode, defectMode);
				impuls_residual = Ddot(2 * N_M + N_P, defectMode, defectMode);
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
						IntoL20FEFunction(modeSolution_i + 2 * N_M, N_P, PressureMode_FeSpace, velocity_space_code_mode, pressure_space_code_mode);

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
					memset(defectMode, 0, (2 * N_M + N_P) * SizeOfDouble);

					SystemMatrixModeAll[subSpaceNum]->GetTNSEResidual(modeSolution_i, defectMode);

					// correction due to L^2_O Pressure space
					if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
						IntoL20Vector2D(defectMode + 2 * N_M, N_P, pressure_space_code_mode);

					residual = Ddot(2 * N_M + N_P, defectMode, defectMode);
					impuls_residual = Ddot(2 * N_M + N_P, defectMode, defectMode);
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
				for (int i = 0; i < N_M; i++)
				{
					solMode1[j * N_M + i] = solModeAll[(2 * j * N_M) + i];
				}
			}
			reorthonormalizeC(solMode1, N_M, subDim);
			for (int j = 0; j < subDim; j++)
			{
				for (int i = 0; i < N_M; i++)
				{
					solModeAll[(2 * j * N_M) + i] = solMode1[j * N_M + i];
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
		fileoutCoeff = generateFileName(coeffBaseName, m, N_Realisations);
		printToTxt(fileoutCoeff, CoeffVector, N_Realisations, subDim, 'C');

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

	for (int sd = 0; sd < subDim; sd++)
	{

		std::string filenameMode = "Mode_NRealisations_" + std::to_string(N_Realisations) + "_ModeN0_" + std::to_string(sd);
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

	CloseFiles();

	return 0;
} // end main
