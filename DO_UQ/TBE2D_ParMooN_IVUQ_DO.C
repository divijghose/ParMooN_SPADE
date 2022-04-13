// =======================================================================
// Purpose:     main program for solving a time-dependent Burger's equation in ParMooN
// Author:      Sashikumaar Ganesan
// History:     Implementation started on 28.11.2020
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <FESpace2D.h>
#include <SystemTBE2D.h>
#include <SystemTBE_Mode2D.h>
#include <SquareStructure2D.h>
#include <Output2D.h>
#include <MainUtilities.h>
#include <TimeDiscRout.h>
#include <TNSE2D_ParamRout.h>

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
#include "../Examples/DO_UQ/burgers_do_test.h" // smooth sol in unit square

// =======================================================================
// main program
// =======================================================================
int main(int argc, char *argv[])
{
	// ======================================================================
	//  declaration of variables
	// ======================================================================
	int i, j, l, m, N_Cells, ORDER, N_U, N_M, N_L, N_TotalDOF, img = 1, N_Modes, N_Realiz;
	int N_Total_MeanDOF, Max_It, NSEType, N_SubSteps, Disctype;

	double *sol, *rhs, *oldrhs, *defect, t1, t2, residual, impuls_residual;
	double limit, AllErrors[7], end_time, oldtau, tau;

	double *solMean, *rhsMean, *old_rhsMean, *old_solMean;
	double *solMode, *rhsMode, *old_rhsMode, *old_solMode;

	TDomain *Domain;
	TDatabase *Database = new TDatabase();
	TFEDatabase2D *FEDatabase = new TFEDatabase2D();
	TCollection *coll, *mortarcoll = NULL;

	TFESpace2D *Velocity_FeSpace, *VelocityMode_FeSpace, *VelocityMean_FeSpace;
	TFESpace2D *fesp[2], *fesp_mode[2], *fesp_mean[2];
	TFESpace2D *Velocity_FeSpace_Mean, *fespMean[2];
	TFESpace2D *Velocity_FeSpace_Mode, *fespMode[2];
	// TFESpace2D *Velocity_FeSpace, *fespMean[2];

	// TFEVectFunct2D *Velocity_Mean, *Velocity_Mode;
	TFEVectFunct2D *Velocity_Mean, *Velocity_Mode, *Velocity;
	// TFEFunction2D *u1_mean, *u2_mean, *fefct[2];

	TFEFunction2D *u1, *u2, *fefct[2];
	TFEFunction2D *u1Mean, *u2Mean, *fefctMean[2];
	TFEFunction2D *u1Mode, *u2Mode, *fefctMode[2];

	TOutput2D *Output, *OutputMean, *OutputMode;

	// TSystemTBE2D *SystemMatrix_Mean;
	TSystemTBE2D *SystemMatrix, *SystemMatrix_Mean, *SystemMatrix_Mode;

	// TSystemTBE_Mode2D *SystemMatrix_Mode;
	// TFEFunction2D *FeFct[2], *FeFct_Mode[4];

	TFEFunction2D *FeFct[2];

	TAuxParam2D *BEaux, *BEaux_error;
	TAuxParam2D *BEaux_mean, *BEaux_error_mean;
	TAuxParam2D *BEaux_mode, *BEaux_error_mode;

	const char vtkdir[] = "VTK";

	char *PsBaseName, *VtkBaseName, *GEO;

	char *VtkBaseNameMean, *VtkBaseNameMode;
	// char UString[] = "u_mean";
	// char UString[] = "usol";
	// char UMString[] = "u_mode";
	// char NameString[] = "UQ";

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
	Disctype = TDatabase::ParamDB->DISCTYPE;

	N_Modes = 2;
	// N_Realiz = TDatabase::ParamDB->P11;

	// VelocityMode = new TFEVectFunct2D*[N_Modes];

	coll = Domain->GetCollection(It_Finest, 0);
	N_Cells = coll->GetN_Cells();
	OutPut("N_Cells : " << N_Cells << endl);

	Velocity_FeSpace = new TFESpace2D(coll, (char *)"Sol", (char *)"FESpace for Solution", BoundCondition, ORDER, NULL);

	VelocityMode_FeSpace = new TFESpace2D(coll, (char*)"Mode", (char*)"Mode", BoundCondition, ORDER, NULL);

	N_U = Velocity_FeSpace->GetN_DegreesOfFreedom();
	N_M = VelocityMode_FeSpace->GetN_DegreesOfFreedom();
	N_Total_MeanDOF = 2*N_U;
	N_TotalDOF = 2*(N_U + N_Modes*N_M);
	// N_TotalDOF = 2 * N_U;
	OutPut("Dof Mean Velocity : "<< setw(10) << 2* N_U << endl);
	OutPut("Dof Mode Velocity : "<< setw(10) << 2* N_M *N_Modes<< endl);
	OutPut("Total Dof all: "<< setw(10) << N_TotalDOF  << endl);
	// OutPut("DOF : " << setw(10) << 2 * N_U << endl);
	exit(0);

	// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	///////////////////////////////////////////////////////////////////////////////////////////////
	////////// -------- REALISATION DATA GENERATION ----------------------------------------- //////
	///////////////////////////////////////////////////////////////////////////////////////////////
	int N_Realisations = TDatabase::ParamDB->REALIZATIONS;
	double LengthScale = TDatabase::ParamDB->LENGTHSCALE;
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

	double *C = new double[N_U * N_U];	// MATRIX
	double *C1 = new double[N_U * N_U]; // MATRIX  - Corelation Matrix
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

			if (TDatabase::ParamDB->stddev_switch == 0)
			{
				double sig_r1 = exp(-1.0 / (1.0 - pow((2 * actual_x - 1), 4))) * exp(-1.0 / (1 - pow((2 * actual_y - 1), 4)));
				double sig_r2 = exp(-1.0 / (1.0 - pow((2 * local_x - 1), 4))) * exp(-1.0 / (1 - pow((2 * local_y - 1), 4)));

				// Co Variance
				C[j * N_U + i] *= sig_r1 * sig_r2 * 5.0;
			}

			else if (TDatabase::ParamDB->stddev_switch == 1)
			{
				double E = TDatabase::ParamDB->stddev_denom;
				double disp = TDatabase::ParamDB->stddev_disp;
				double power = TDatabase::ParamDB->stddev_power;
				double sig_r1 = exp(-pow((2 * actual_x - 1 - disp), power) / (E)) / (2 * 3.14159265359 * sqrt(E)) * exp(-pow((2 * actual_x - 1 - disp), power) / (E)) / (2 * 3.14159265359 * sqrt(E));
				double sig_r2 = exp(-pow((2 * local_x - 1 - disp), power) / (E)) / (2 * 3.14159265359 * sqrt(E)) * exp(-pow((2 * local_y - 1 - disp), power) / (E)) / (2 * 3.14159265359 * sqrt(E));
				// Co Variance
				C[j * N_U + i] *= sig_r1 * sig_r2;
			}

			else
			{
				cout << "Error " << endl;
				exit(0);
			}

			norm += C[j * N + i] * C[j * N + i];
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

	info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', m1, n, C, lda,
						  S, U, ldu, Vt, ldvt, superb);

	cout << endl
		 << endl;

	if (info > 0)
	{
		printf("The algorithm computing SVD failed to converge.\n");
		exit(1);
	}
	cout << " SVD COMPUTED " << endl;
	int energyVal = 0;

	double sumSingularVal = 0;
	for (int i = 0; i < N_U; i++)
		sumSingularVal += S[i];
	double val = 0;
	for (energyVal = 0; energyVal < N_U; energyVal++)
	{
		val += S[energyVal];
		if (val / sumSingularVal > EigenPercent)
			break;
	}

	cout << " MODES : " << energyVal + 1 << endl;

	int modDim = energyVal + 1;

	double *Ut = new double[N_U * modDim]();
	double *Z = new double[N_Realisations * modDim]();

	double *RealizationVector = new double[N_U * N_Realisations]();
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
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_U, N_Realisations, modDim, 1.0, Ut, modDim, Z, N_Realisations, 0.0, RealizationVector, N_Realisations);
	cout << " MULT DONE " << endl;
	// printMatrix(RealizationVector, N_DOF,N_Realisations);

	// mkl_dimatcopy('R','T', N_DOF,N_Realisations,1.0,RealizationVector,N_DOF,N_Realisations);
	cout << " COPY DONE " << endl;

	cout << " REALISATIONS COMPUTED " << endl;

	//////////////////////////////////End of Realization/////////////////////////////////////////

	srand(time(NULL));
	int N_samples = 100;
	int *indexArray = new int[N_samples];
	for (int i = 0; i < N_samples; i++)
		indexArray[i] = rand() % N_U;

	double *RealizationVectorTemp = new double[N_U * N_Realisations]();
	for (int i = 0; i < N_U; i++)
	{
		for (int j = 0; j < N_Realisations; j++)
		{
			RealizationVectorTemp[mappingArray[i] * N_Realisations + j] = RealizationVector[j + N_Realisations * i];
		}
	}

	memcpy(RealizationVector, RealizationVectorTemp, N_U * N_Realisations * SizeOfDouble);

	/////////////////////////////////////// -------- END OF REALISATION DATA SETS ------------ ////////////////////////////////////////////////////////////////

	////////////////////////////////////// -------- START OF DO INITIALIZATION ------------ ////////////////////////////////////////////////////////////////

	double *MeanVector = new double[N_U * 1](); // overline{C}_{dof} = \sum_{i=1}^{N_Realisations}(C^{i}_{dof})/N_Realisations
	for (int i = 0; i < N_U; ++i)
	{
		for (int j = 0; j < N_Realisations; ++j)
		{
			MeanVector[i] += (RealizationVector[i * N_Realisations + j] / N_Realisations);
		}
	}

	double *PerturbationVector = new double[N_U * N_Realisations](); // \hat{C}^{i}_{dof} = C^{i}_{dof} - \overline{C}_{dof}
	for (int i = 0; i < N_U; ++i)
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
	cout << " DO SVD COMPUTED " << endl;

	//////////////////////////////////////////// DO - SVD End///////////////////////////////

	///////DO - Subspace dimension calculation /////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////
	double SVPercent = TDatabase::ParamDB->SVPERCENT;
	int s = 0;
	double valDO = 0.0;
	double sumSingularValDO = 0.0;
	for (int i = 0; i < minDim; i++)
	{
		sumSingularValDO += Sg[i];
	}
	while (valDO / sumSingularValDO < SVPercent)
	{
		valDO += Sg[s];
		s++;
	}

	cout << " SUBSPACE DIMENSION : " << s + 1 << endl;

	int subDim = s + 1;
	subDim=1;
	////////Subspace dimension calculated//////////////////
	TDatabase::ParamDB->N_Subspace_Dim = subDim;

	/////Projection Matrix///////////
	////////////////////////////////
	cout << " Min DIMENSION : " << minDim << endl;
	double *ProjectionVector = new double[N_Realisations * minDim]();
	cout << "PROJ VECTOR MULT START " << endl;
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N_Realisations, minDim, N_U, 1.0, PerturbationVector, N_Realisations, L, minDim, 0.0, ProjectionVector, minDim);
	cout << "PROJ VECTOR MULT DONE " << endl;

	/// Initialize Coefficient Matrix - First subDim columns of Projection Matrix ////////////////////
	double *CoeffVector = new double[N_Realisations * subDim]();
	// memcpy(CoeffVector, ProjectionVector, N_Realisations*subDim*SizeOfDouble); //For ColMajor storage -wrong!!
	// for (int i=0;i<N_Realisations;i++){
	// 	for (int j=0;j<subDim;j++){
	// 		CoeffVector[i*subDim+j] = ProjectionVector[i*minDim+j];
	// 	}
	// }
	for (int i = 0; i < N_Realisations; i++)
	{
		for (int j = 0; j < subDim; j++)
		{
			CoeffVector[j * N_Realisations + i] = ProjectionVector[i * minDim + j]; // CoeffVector in Col Major form
		}
	}

	////////////Initialize Mode Vector - First subDim columns of Left Singular Vector//////////////////
	double *ModeVector = new double[N_U * subDim]();
	// memcpy(ModeVector, L, N_U*subDim*SizeOfDouble);//For ColMajor storage
	// for (int i=0;i<N_U;i++){
	// 	for (int j=0;j<subDim;j++){
	// 		ModeVector[i*subDim+j] = ProjectionVector[i*minDim+j];
	// 	}
	// }
	for (int i = 0; i < N_U; i++)
	{
		for (int j = 0; j < subDim; j++)
		{
			ModeVector[j * N_U + i] = L[i * minDim + j]; // ModeVector in Col Major form
		}
	}

	////////////////////////////////////////////DO - Initialization Ends//////////////////////////////////////
	///////================================================================================//////////////////

	// xxxxxxxxxxxxxxxxxxxxxxxxxxx

	//

	//======================================================================
	// construct all finite element functions
	//======================================================================
	sol = new double[N_TotalDOF];
	rhs = new double[N_TotalDOF];
	oldrhs = new double[N_TotalDOF];
	defect = new double[N_TotalDOF];

	solMean = new double[N_TotalDOF];
	old_solMean = new double[N_TotalDOF];
	rhsMean = new double[N_TotalDOF];
	old_rhsMean = new double[N_TotalDOF];

	solMode = new double[N_TotalDOF * subDim]();
	rhsMode = new double[N_TotalDOF * subDim]();
	old_rhsMode = new double[N_TotalDOF]();

	memset(sol, 0, N_TotalDOF * SizeOfDouble);
	memset(rhs, 0, N_TotalDOF * SizeOfDouble);
	memset(defect, 0, N_TotalDOF * SizeOfDouble);

	memset(solMean, 0, N_TotalDOF * SizeOfDouble);
	memset(old_solMean, 0, N_TotalDOF * SizeOfDouble);
	memset(rhsMean, 0, N_TotalDOF * SizeOfDouble);
	memset(old_rhsMean, 0, N_TotalDOF * SizeOfDouble);

	memset(solMode, 0, N_TotalDOF * subDim * SizeOfDouble);
	memset(rhsMode, 0, N_TotalDOF * subDim * SizeOfDouble);
	memset(old_rhsMode, 0, N_TotalDOF * SizeOfDouble);

	Velocity_Mean = new TFEVectFunct2D(Velocity_FeSpace, (char *)"U_Mean", (char *)"Mean Component", solMean, N_U, 2); // check length
	Velocity_Mode = new TFEVectFunct2D(Velocity_FeSpace, (char *)"U_Mode", (char *)"Mode", solMode, N_TotalDOF, subDim);

	u1Mean = Velocity_Mean->GetComponent(0);
	u2Mean = Velocity_Mean->GetComponent(1);

	// u1Mean->Interpolate(InitialU1Mean);
	// u2Mean->Interpolate(InitialU2Mean);

	for (i = 0; i < N_U; i++)
	{
		// sol[mappingArray[i]] = RealizationVector[RealNo+N_Realisations*i]/10;
		// sol[N_U+mappingArray[i]]=RealizationVector[RealNo+N_Realisations*i]/10;
		solMean[i] = MeanVector[i];
		solMean[N_U + i] = MeanVector[i]; // no second component
	}

	//======================================================================
	// /DO - SystemMatrix construction and solution
	//======================================================================
	// Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
	// Solver: AMG_SOLVE (or) GMG  (or) DIRECT
	SystemMatrix_Mean = new TSystemTBE2D(Velocity_FeSpace, Velocity_Mean, solMean, rhsMean, GALERKIN, DIRECT);

	// SystemMatrix_Mode = new TSystemTBE2D(Velocity_FeSpace, Velocity_Mode, solMode, rhsMode, GALERKIN, DIRECT);

	fefctMean[0] = u1Mean;
	fefctMean[1] = u2Mean;
	fespMean[0] = Velocity_FeSpace;
	fespMode[0] = Velocity_FeSpace;


	BEaux_mean = new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2, TimeNSN_ParamFct2,
							TimeNSN_FEValues2, fespMean, FeFct, TimeNSFct2, TimeNSFEFctIndex2,
							TimeNSFEMultiIndex2, TimeNSN_Params2, TimeNSBeginParam2);

	// aux for calculating the error
	if (TDatabase::ParamDB->MEASURE_ERRORS)
	{
		BEaux_error_mean = new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
									  TimeNSN_ParamFct2,
									  TimeNSN_FEValues2,
									  fespMean, FeFct,
									  TimeNSFct2,
									  TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
									  TimeNSN_Params2, TimeNSBeginParam2);
	}

	BEaux_mode = new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2, TimeNSN_ParamFct2,
							TimeNSN_FEValues2, fespMean, FeFct, TimeNSFct2, TimeNSFEFctIndex2,
							TimeNSFEMultiIndex2, TimeNSN_Params2, TimeNSBeginParam2);

	// aux for calculating the error
	if (TDatabase::ParamDB->MEASURE_ERRORS)
	{
		BEaux_error_mode = new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
									  TimeNSN_ParamFct2,
									  TimeNSN_FEValues2,
									  fespMode, FeFct,
									  TimeNSFct2,
									  TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
									  TimeNSN_Params2, TimeNSBeginParam2);
	}

	// initilize the system matrix with the functions defined in Example file
	SystemMatrix_Mean->Init(DO_Mean_Equation_Coefficients, BoundCondition, U1BoundValue, U2BoundValue, BEaux_mean, BEaux_error_mean);
    SystemMatrix_Mean->Assemble(solMean, rhsMean);


	// SystemMatrix_Mode->Init(DO_Mode_Equation_Coefficients, BoundCondition, U1BoundValue, U2BoundValue, BEaux_mode, BEaux_error_mode);

	// -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
	//------------------------------------------ MEAN EQUATION SETUP -----------------------------------------------------//
	// -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//

	// Aux Setup for the RHS -- There is no Aux for the Mean equation, So set the values as NULL
	TAuxParam2D *aux;
	// fespMean[0] = Velocity_FeSpace;
	aux = new TAuxParam2D(1, 0, 0, 0, fespMean, NULL, NULL, NULL, NULL, 0, NULL);

	/* -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-
	--------------------------------------[[[ END  ]]] MEAN EQUATION SETUP -----------------------------------------------------
	 -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-*/
	//  -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
	//------------------------------------------ CO EFFICIENT EQUATION SETUP -----------------------------------------------------//
	// -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//

	TFEVectFunct2D *FeVector_Coefficient = new TFEVectFunct2D(Velocity_FeSpace, (char *)"Coeff", (char *)"Coefficients", CoeffVector, N_Realisations, subDim);

	// -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
	//------------------------------------------ CO EFFICIENT EQUATION SETUP END -------------------------------------------------------//
	// -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//

	// -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
	//-------------------------------------- MODE EQUATION SETUP -----------------------------------------------------//
	// -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
	// int N_terms_Mode = 3; // Number of Shape function derivatives required ( in this case 3, N, NX, NY )
	// MultiIndex2D Derivatives_MatrixARhs_Mode[3] = {D00, D10, D01};
	// int SpacesNumbers_MatrixARhs_Mode[3] = {0, 0, 0};
	// int N_Matrices_MatrixARhs_Mode = 1;
	// int RowSpace_MatrixARhs_Mode[1] = {0};
	// int ColumnSpace_MatrixARhs_Mode[1] = {0};
	// int N_Rhs_MatrixARhs_Mode = 1;
	// int RhsSpace_MatrixARhs_Mode[1] = {0};
	// SystemMatrix_Mode->Init_WithDiscreteform(DO_Mode_Equation_Coefficients, BoundCondition, BoundValue, "DO_LINEAR_Mode", "DO_LINEAR_Mode",
	//                                          N_terms_Mode, Derivatives_MatrixARhs_Mode, SpacesNumbers_MatrixARhs_Mode,
	//                                          N_Matrices_MatrixARhs_Mode, N_Rhs_MatrixARhs_Mode,
	//                                          RowSpace_MatrixARhs_Mode, ColumnSpace_MatrixARhs_Mode, RhsSpace_MatrixARhs_Mode,
	//                                          DO_Mode_Equation_Assembly, DO_Mode_Equation_Coefficients,
	//                                          NULL);

	for (int j = 0; j < subDim; j++)
	{
		for (int i = 0; i < N_U; i++)
		{
			// solMode[j*N_DOF+mappingArray[i]] = ModeVector[j*N_DOF+i];
			solMode[j * N_TotalDOF + i] = ModeVector[j * N_U + i];
			solMode[j * N_TotalDOF + N_U + i] = ModeVector[j * N_U + i];
		}
	}

	CalcCovarianceMatx(CoeffVector);
	CalcCoskewnessMatx(CoeffVector);
	InvertCov();

	// 	// double* ModeVector_OldRHS = new double[N_DOF]();

	// // Set up a FE VECT FUNCTION TO STORE ALL THE Components of CTilde

	//     // TFEVectFunct2D* linModesFeVectFunct =

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

	// ;
	//     // Set Up Aux  Param for the given MOde
	//     // The mode equation needs the entire values of the C_tilde matrix array
	//     // TAuxParam2D* aux_RHS_DO = new TAuxParam2D (	TimeLinear_FESpaces_DO, <FE VECT FUNCTION>, TimeLinear_ParamFct_DO,
	// 	// 											TimeLinear_FEValues_DO,
	// 	// 											fesp_RHS,
	// 	// 											TimeNSFct_DO,
	// 	// 											TimeNSFEMultiIndex_DO,
	// 	// 											TimeLinear_Params_DO, TimeNSBeginParam_DO);

	TAuxParam2D *aux_RHS_DO = new TAuxParam2D(1, 0, 0, 0, fespMean, NULL, NULL, NULL, NULL, 0, NULL);
	// -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
	//--------------------------------------[[[ END  ]]] MODE EQUATION SETUP -----------------------------------------------------//
	// -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//

	// // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	//======================================================================
	// produce outout at t=0
	//======================================================================
	// VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
	std::string strMean = std::to_string(N_Realisations);
	std::string filenameMean = "Mean_NRealisations_" + std::to_string(N_Realisations);
	VtkBaseNameMean = const_cast<char *>(filenameMean.c_str());

	std::string strMode = std::to_string(N_Realisations);
	std::string filenameMode = "Mode_NRealisations_" + std::to_string(N_Realisations);
	VtkBaseNameMode = const_cast<char *>(filenameMode.c_str());

	std::string fileoutCoeff;
	std::string fileoutMean;
	std::string fileoutMode;
	std::string fileoutMC;

	OutputMean = new TOutput2D(2, 2, 1, 1, Domain);
	OutputMean->AddFEVectFunct(Velocity_Mean);

	OutputMode = new TOutput2D(2, 2, 1, 1, Domain);
	OutputMode->AddFEVectFunct(Velocity_Mode);

	int meanimg = 0;
	int modeimg = 0;

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

	if (TDatabase::ParamDB->WRITE_VTK)
	{
		os.seekp(std::ios::beg);
		if (modeimg < 10)
			os << "VTK/" << VtkBaseNameMode << ".0000" << modeimg << ".vtk" << ends;
		else if (modeimg < 100)
			os << "VTK/" << VtkBaseNameMode << ".000" << modeimg << ".vtk" << ends;
		else if (modeimg < 1000)
			os << "VTK/" << VtkBaseNameMode << ".00" << modeimg << ".vtk" << ends;
		else if (modeimg < 10000)
			os << "VTK/" << VtkBaseNameMode << ".0" << modeimg << ".vtk" << ends;
		else
			os << "VTK/" << VtkBaseNameMode << "." << modeimg << ".vtk" << ends;
		OutputMode->WriteVtk(os.str().c_str());
		modeimg++;
	}

	// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	// coll->GetHminHmax(&hmin, &hmax);
	// OutPut("h_min : " << hmin << " h_max : " << hmax << endl);

	//======================================================================
	// time disc loop
	//======================================================================
	// parameters for time stepping scheme
	m = 0;
	fileoutMean = "Mean/Mean_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m);
	std::ofstream fileMean;
	fileMean.open(fileoutMean);

	for (int i = 0; i < N_TotalDOF; i++)
	{

		fileMean << MeanVector[i] << ",";
		fileMean << endl;
	}

	fileMean.close();

	fileoutMode = "Modes/Mode_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m);
	std::ofstream fileMode;
	fileMode.open(fileoutMode);

	for (int i = 0; i < N_TotalDOF; i++)
	{
		for (int j = 0; j < subDim; j++)
		{
			fileMode << solMode[j * subDim + i];
			if (j != subDim - 1)
				fileMode << ",";
		}
		fileMode << endl;
	}

	fileMode.close();

	fileoutCoeff = "Coefficients/Coeff_NRealisations_" + std::to_string(N_Realisations) + "_t" + std::to_string(m);
	std::ofstream fileCoeff;
	fileCoeff.open(fileoutCoeff);

	for (int i = 0; i < N_Realisations; i++)
	{
		for (int j = 0; j < subDim; j++)
		{
			fileCoeff << CoeffVector[j * subDim + i];
			if (j != subDim - 1)
				fileCoeff << ",";
		}
		fileCoeff << endl;
	}

	fileCoeff.close();
	N_SubSteps = GetN_SubSteps();
	oldtau = 1.;
	end_time = TDatabase::TimeDB->ENDTIME;
	limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
	Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
	memset(AllErrors, 0, 7. * SizeOfDouble);

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

			// OutPut(endl << "CURRENT TIME: ");
			// OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

			// copy sol, rhs to olssol, oldrhs
			// memcpy(old_rhsMean, rhsMean, N_TotalDOF * SizeOfDouble);
			// memcpy(old_solMean, solMean, N_TotalDOF * SizeOfDouble);

			// DO_Mean_RHS(Velocity_FeSpace, Velocity_Mode,subDim, rhsMean,N_U);

			// assemble only rhs, nonlinear matrix for NSE will be assemble in fixed point iteration
			// not needed if rhs is not time-dependent
			if (m != 1)
			{
				SystemMatrix_Mean->AssembleA();

				//  SystemMatrix->AssembleA();
			}
			else
			{
				//  SystemMatrix_Mean->Assemble(solMean, rhsMean);
				SystemMatrix_Mean->Assemble(solMean, rhsMean);
			}

			cout << "Comes Here" << endl;
			exit(0);

			// scale B matices and assemble NSE-rhs based on the \theta time stepping scheme
			//  SystemMatrix_Mean->AssembleSystMat(oldrhs, rhs, sol);
			SystemMatrix_Mean->AssembleSystMat(old_rhsMean, rhsMean, solMean);
			oldtau = tau;

			// calculate the residual
			// SystemMatrix_Mean->GetTBEResidual(sol, defect);
			SystemMatrix_Mean->GetTBEResidual(solMean, defect);

			residual = Ddot(N_TotalDOF, defect, defect);
			// OutPut("Nonlinear iteration step   0");
			// OutPut(setw(14) << sqrt(residual) << endl);

			//======================================================================
			// Solve the system
			// Nonlinear iteration of fixed point type
			//======================================================================
			for (j = 1; j <= Max_It; j++)
			{
				// Solve the NSE system
				//  SystemMatrix_Mean->Solve(sol);
				SystemMatrix_Mean->Solve(solMean);

				// restore the mass matrix for the next nonlinear iteration
				//  SystemMatrix_Mean->RestoreMassMat();
				SystemMatrix_Mean->RestoreMassMat();

				// assemble the system matrix with given aux, sol and rhs
				//  SystemMatrix_Mean->AssembleANonLinear(sol, rhs);
				SystemMatrix_Mean->AssembleANonLinear(solMean, rhsMean);

				// assemble system mat, S = M + dt\theta_1*A
				//  SystemMatrix_Mean->AssembleSystMatNonLinear();
				SystemMatrix_Mean->AssembleSystMatNonLinear();

				// get the residual
				memset(defect, 0, N_TotalDOF * SizeOfDouble);
				//  SystemMatrix_Mean->GetTBEResidual(sol, defect);
				SystemMatrix_Mean->GetTBEResidual(solMean, defect);

				residual = Ddot(N_TotalDOF, defect, defect);
				//  OutPut("nonlinear iteration step " << setw(3) << j);
				//  OutPut(setw(14) << sqrt(residual) << endl);

				if (sqrt(residual) <= limit)
					break;

			} // for(j=1;j<=Max_It;j++)

			// restore the mass matrix for the next time step
			// SystemMatrix_Mean->RestoreMassMat();
			SystemMatrix_Mean->RestoreMassMat();

		} // for(l=0;l<N_SubSteps;
		  //======================================================================
		  // measure errors to known solution
		  //======================================================================
		if (TDatabase::ParamDB->MEASURE_ERRORS)
		{
			// SystemMatrix_Mean->MeasureErrors(ExactU1, ExactU2, AllErrors);
			SystemMatrix_Mean->MeasureErrors(ExactU1, ExactU2, AllErrors);

			//  OutPut("L2(u): " <<   AllErrors[0] << endl);
			//  OutPut("H1-semi(u): " <<  AllErrors[1] << endl);
			//  OutPut("L2(p): " <<  AllErrors[2] << endl);
			//  OutPut("H1-semi(p): " <<  AllErrors[3]   << endl);
			//  OutPut(AllErrors[4] <<  " l_infty(L2(u)) " <<AllErrors[5] << endl);
			//  OutPut(TDatabase::TimeDB->CURRENTTIME << " L2(0,t,L2)(u) : " <<   sqrt(AllErrors[6]) << endl);
		} // if(TDatabase::ParamDB->MEASURE_ERRORS)

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

		for (int subSpaceNum = 0; subSpaceNum < subDim; subSpaceNum++)
		{

			DO_CoEfficient(Velocity_FeSpace, Velocity_Mode, FeVector_Coefficient, Velocity_Mean, subDim, subSpaceNum, N_Realisations);
		}

		CalcCovarianceMatx(CoeffVector);
		CalcCoskewnessMatx(CoeffVector);
		InvertCov();

		// xxxxxxxxxxxxxxxxxxxxxxxxxxxx
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

			for (int subSpaceNum = 0; subSpaceNum < subDim; subSpaceNum++)
			{																 // subspace loop
				double *modeSolution_i = solMode + subSpaceNum * N_TotalDOF; // This works for column major
				double *modeSolution_rhs = rhsMode + subSpaceNum * N_TotalDOF;

				// copy sol, rhs to olssol, oldrhs
				memcpy(old_rhsMode, modeSolution_rhs, N_TotalDOF * SizeOfDouble);

				// Assemble rhs
				DO_Mode_RHS(Velocity_FeSpace, Velocity_Mean, Velocity_Mode, subDim, modeSolution_rhs, subSpaceNum);

				SystemMatrix_Mode->Assemble(modeSolution_i, modeSolution_rhs);
				//   }

				// scale B matices and assemble NSE-rhs based on the \theta time stepping scheme
				//  SystemMatrix_Mean->AssembleSystMat(oldrhs, rhs, sol);
				SystemMatrix_Mode->AssembleSystMat(old_rhsMode, modeSolution_rhs, modeSolution_i);
				oldtau = tau;

				// calculate the residual
				// SystemMatrix_Mean->GetTBEResidual(sol, defect);
				SystemMatrix_Mode->GetTBEResidual(modeSolution_i, defect);

				residual = Ddot(N_TotalDOF, defect, defect);
				// OutPut("Nonlinear iteration step   0");
				// OutPut(setw(14) << sqrt(residual) << endl);

				//======================================================================
				// Solve the system
				// Nonlinear iteration of fixed point type
				//======================================================================
				for (j = 1; j <= Max_It; j++)
				{
					// Solve the NSE system
					//  SystemMatrix_Mean->Solve(sol);
					SystemMatrix_Mode->Solve(modeSolution_i);

					// restore the mass matrix for the next nonlinear iteration
					//  SystemMatrix_Mean->RestoreMassMat();
					SystemMatrix_Mode->RestoreMassMat();

					// assemble the system matrix with given aux, sol and rhs
					//  SystemMatrix_Mean->AssembleANonLinear(sol, rhs);
					SystemMatrix_Mode->AssembleANonLinear(modeSolution_i, modeSolution_rhs);

					// assemble system mat, S = M + dt\theta_1*A
					//  SystemMatrix_Mean->AssembleSystMatNonLinear();
					SystemMatrix_Mode->AssembleSystMatNonLinear();

					// get the residual
					memset(defect, 0, N_TotalDOF * SizeOfDouble);
					//  SystemMatrix_Mean->GetTBEResidual(sol, defect);
					SystemMatrix_Mode->GetTBEResidual(modeSolution_i, defect);

					residual = Ddot(N_TotalDOF, defect, defect);
					//  OutPut("nonlinear iteration step " << setw(3) << j);
					//  OutPut(setw(14) << sqrt(residual) << endl);

					if (sqrt(residual) <= limit)
						break;

				} // for(j=1;j<=Max_It;j++)
				SystemMatrix_Mode->RestoreMassMat();
			} // subspace loop end

			// restore the mass matrix for the next time step
			// SystemMatrix_Mean->RestoreMassMat();
			// SystemMatrix_Mean->RestoreMassMat();

		} // l substep time loop
		  //======================================================================
		// produce outout
		//======================================================================
		if (m == 1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
		{

			if (TDatabase::ParamDB->WRITE_VTK)
			{
				os.seekp(std::ios::beg);
				if (modeimg < 10)
					os << "VTK/" << VtkBaseNameMode << ".0000" << modeimg << ".vtk" << ends;
				else if (modeimg < 100)
					os << "VTK/" << VtkBaseNameMode << ".000" << modeimg << ".vtk" << ends;
				else if (modeimg < 1000)
					os << "VTK/" << VtkBaseNameMode << ".00" << modeimg << ".vtk" << ends;
				else if (modeimg < 10000)
					os << "VTK/" << VtkBaseNameMode << ".0" << modeimg << ".vtk" << ends;
				else
					os << "VTK/" << VtkBaseNameMode << "." << modeimg << ".vtk" << ends;
				OutputMode->WriteVtk(os.str().c_str());
				modeimg++;
			}
		}

		// xxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	} // while(TDatabase::TimeDB->CURRENTTIME< e

	//======================================================================
	// produce final outout
	//======================================================================
	if (TDatabase::ParamDB->WRITE_VTK)
	{
		os.seekp(std::ios::beg);
		if (img < 10)
			os << "VTK/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
		else if (img < 100)
			os << "VTK/" << VtkBaseName << ".000" << img << ".vtk" << ends;
		else if (img < 1000)
			os << "VTK/" << VtkBaseName << ".00" << img << ".vtk" << ends;
		else if (img < 10000)
			os << "VTK/" << VtkBaseName << ".0" << img << ".vtk" << ends;
		else
			os << "VTK/" << VtkBaseName << "." << img << ".vtk" << ends;
		Output->WriteVtk(os.str().c_str());
		img++;
	}
	TDatabase::TimeDB->CURRENTTIME = 0;

	CloseFiles();

	return 0;
} // end main
