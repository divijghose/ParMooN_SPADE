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
#include "../Examples/DO_UQ/burgers_do_validation.h" // smooth sol in unit square

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
	// TFEFunction2D *u1_mean, *u2_mean, *fefct[2];

	TFEFunction2D *u1, *u2, *fefct[2];
	TFEFunction2D *u1Mean, *u2Mean, *fefctMean[2];
	TFEFunction2D *u1Mode, *u2Mode, *fefctMode[2];

	// Triple pointer for sending in a FEVect function for each of the Aux
	TFEFunction2D *fefctModeAll[50][2];

	TOutput2D *Output, *OutputMean, *OutputMode;

	// TSystemTBE2D *SystemMatrix_Mean;

	// TSystemTBE_Mode2D *SystemMatrix_Mode;
	// TFEFunction2D *FeFct[2], *FeFct_Mode[4];

	TFEFunction2D *FeFct[2];

	TAuxParam2D *BEaux, *BEaux_error;
	TAuxParam2D *BEaux_mean, *BEaux_error_mean;
	TAuxParam2D *BEaux_mode, *BEaux_error_mode;

	TAuxParam2D **BEaux_modeAll;

	char *PsBaseName, *VtkBaseName, *GEO;

	char *VtkBaseNameMean, *VtkBaseNameMode;
	// char UString[] = "u_mean";
	// char UString[] = "usol";
	// char UMString[] = "u_mode";
	// char NameString[] = "UQ";

	std::ostringstream os;
	os << " ";

	const char vtkdir[] = "VTK";
	const char modedir[] = "Modes";
	const char meandir[] = "Mean";
	const char coeffdir[] = "Coefficients";
	const char mcdir[] = "MonteCarlo";
	const char endir[] = "Energy_Data";

	mkdir(vtkdir, 0777);
	mkdir(meandir, 0777);
	mkdir(modedir, 0777);
	mkdir(coeffdir, 0777);
	mkdir(mcdir, 0777);
	mkdir(endir, 0777);

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

	coll = Domain->GetCollection(It_Finest, 0);
	N_Cells = coll->GetN_Cells();
	OutPut("Number of Cells : " << N_Cells << endl);

	Velocity_FeSpace = new TFESpace2D(coll, (char *)"Sol", (char *)"FESpace for Solution", BoundCondition, ORDER, NULL);
	N_U = Velocity_FeSpace->GetN_DegreesOfFreedom();

	OutPut("DOF for general solution" << setw(10) << 2 * N_U << endl);

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
					C[i * N_U + j] *= sig_r1 * sig_r2;
				}

				else if (TDatabase::ParamDB->stddev_switch == 2)
				{
					double amplitude = TDatabase::ParamDB->stddev_power;
					double sig_r1 = amplitude * sin(-1.0 * Pi * (2 * actual_x - 2)) * sin(-1.0 * Pi * (2 * actual_y - 2));
					double sig_r2 = amplitude * sin(-1.0 * Pi * (2 * local_x - 2)) * sin(-1.0 * Pi * (2 * local_y - 2));
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

		std::string fileOutCorrelation = "Correlation_" + std::to_string(N_U) + "_NR_" + std::to_string(N_Realisations) + ".txt";
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

		std::string fileOutCorrS = "Corr_S_" + std::to_string(N_U) + "_NR_" + std::to_string(N_Realisations) + ".txt";
		printToTxt(fileOutCorrS, S, N_U, 1, 'C');

		std::string fileOutCorrU = "Corr_U_" + std::to_string(N_U) + "_NR_" + std::to_string(N_Realisations) + ".txt";
		printToTxt(fileOutCorrU, U, N_U, N_U, 'R');

		std::string fileOutCorrVt = "Corr_Vt_" + std::to_string(N_U) + "_NR_" + std::to_string(N_Realisations) + ".txt";
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

	////////////////////////////////////// -------- START OF DO INITIALIZATION ------------ ////////////////////////////////////////////////////////////////

	double *MeanVector = new double[N_U * 1](); // overline{C}_{dof} = \sum_{i=1}^{N_Realisations}(C^{i}_{dof})/N_Realisations
	for (int i = 0; i < N_U; ++i)
	{
		for (int j = 0; j < N_Realisations; ++j)
		{
			MeanVector[i] += (RealizationVector[i * N_Realisations + j] / N_Realisations);
		}
	}
	printToTxt("Mean.txt", MeanVector, N_U, 1, 'C');

	double *PerturbationVector = new double[N_U * N_Realisations](); // \hat{C}^{i}_{dof} = C^{i}_{dof} - \overline{C}_{dof}
	for (int i = 0; i < N_U; ++i)
	{
		for (int j = 0; j < N_Realisations; ++j)
		{
			PerturbationVector[i * N_Realisations + j] = RealizationVector[i * N_Realisations + j] - MeanVector[i];
		}
	}
	printToTxt("PerturbationVector.txt", PerturbationVector, N_U, N_Realisations, 'R');

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

	/////Projection Matrix///////////
	////////////////////////////////
	double *ProjectionVector = new double[N_Realisations * minDim]();

	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N_Realisations, minDim, N_U, 1.0, PerturbationVector, N_Realisations, L, minDim, 0.0, ProjectionVector, minDim);

	/// Initialize Coefficient Matrix - First subDim columns of Projection Matrix ////////////////////
	double *CoeffVector = new double[N_Realisations * subDim]();

	for (int i = 0; i < N_Realisations; i++)
	{
		for (int j = 0; j < subDim; j++)
		{
			CoeffVector[j * N_Realisations + i] = ProjectionVector[i * minDim + j]; // CoeffVector in Col Major form
		}
	}
	printToTxt("Coeff.txt", CoeffVector, N_Realisations, subDim, 'C');

	////////////Initialize Mode Vector - First subDim columns of Left Singular Vector//////////////////

	double *ModeVector = new double[N_U * subDim]();

	for (int i = 0; i < N_U; i++)
	{
		for (int j = 0; j < subDim; j++)
		{
			ModeVector[j * N_U + i] = L[i * minDim + j]; // ModeVector in Col Major form
		}
	}
	printToTxt("Mode.txt", ModeVector, N_U, subDim, 'C');
	const char ip[] = "IPMatrices";
	mkdir(ip, 0777);

	double *IPMatxMode = new double[subDim * subDim]();
	double *IPMatxMean = new double[1 * 1]();

	calcIPMatx(IPMatxMode, ModeVector, N_U, subDim, 'C');
	printToTxt("IPMatxMode_Init.txt", IPMatxMode, subDim, subDim, 'R');

	calcIPMatx(IPMatxMean, MeanVector, N_U, 1, 'C');
	printToTxt("IPMatxMean_Init.txt", IPMatxMean, 1, 1, 'R');

	m = 0;

	// double *ModeVector = new double[N_U * subDim]();

	delete[] PerturbationVector;
	delete[] PerturbationVectorCopy;
	delete[] Sg;
	delete[] L;
	delete[] Rt;
	delete[] ProjectionVector;
	////////////////////////////////////////////DO - Initialization Ends//////////////////////////////////////
	///////================================================================================//////////////////

	VelocityMean_FeSpace = new TFESpace2D(coll, (char *)"Mean", (char *)"FE Space for Mean Solution", BoundCondition, ORDER, NULL);

	VelocityMode_FeSpace = new TFESpace2D(coll, (char *)"Mode", (char *)"FE Space for Mode Solution", BoundCondition, ORDER, NULL);

	N_M = VelocityMode_FeSpace->GetN_DegreesOfFreedom();
	N_Total_MeanDOF = 2 * N_U;
	int N_Total_ModeDOF = 2 * N_M * subDim;
	N_TotalDOF = N_Total_MeanDOF + N_Total_ModeDOF;
	OutPut("Dof Mean Velocity : " << setw(10) << N_Total_MeanDOF << endl);
	OutPut("Dof Mode Velocity : " << setw(10) << N_Total_ModeDOF << endl);
	OutPut("Total DOF Mean+Mode : " << setw(10) << N_TotalDOF << endl);

	//======================================================================
	// construct all finite element functions
	//======================================================================

	double *solMean, *rhsMean, *old_rhsMean, *old_solMean;
	solMean = new double[N_Total_MeanDOF]();
	old_solMean = new double[N_Total_MeanDOF]();
	rhsMean = new double[N_Total_MeanDOF]();
	old_rhsMean = new double[N_Total_MeanDOF]();

	defect = new double[N_Total_MeanDOF]();

	double *solMode, *rhsMode, *old_rhsMode, *old_solMode;
	solMode = new double[N_Total_ModeDOF]();
	old_solMode = new double[N_Total_ModeDOF]();
	rhsMode = new double[N_Total_ModeDOF]();
	old_rhsMode = new double[2 * N_M]();

	double *solModeAll, *rhsModeAll;
	solModeAll = new double[N_Total_ModeDOF]();
	rhsModeAll = new double[N_Total_ModeDOF]();
	for (i = 0; i < N_U; i++)
	{

		solMean[i] = MeanVector[i];
		solMean[N_U + i] = 0;
	}

	for (int j = 0; j < subDim; j++)
	{
		for (int i = 0; i < N_M; i++)
		{
			solMode[(2 * j * N_M) + i] = ModeVector[j * N_U + i];
			solMode[(2 * j * N_M) + N_M + i] = 0;
		}
	}

	for (int j = 0; j < subDim; j++)
	{
		for (int i = 0; i < N_M; i++)
		{
			solModeAll[(2 * j * N_M) + i] = ModeVector[j * N_U + i];
			solModeAll[(2 * j * N_M) + N_M + i] = 0;
		}
	}

	TFEVectFunct2D *Velocity_Mean, *Velocity_Mode, *Velocity;
	Velocity_Mean = new TFEVectFunct2D(VelocityMean_FeSpace, (char *)"U_Mean", (char *)"Mean Component", solMean, N_U, 2); // check length

	Velocity_Mode = new TFEVectFunct2D(VelocityMode_FeSpace, (char *)"U_Mode", (char *)"Mode Component", solMode, N_M, 2 * subDim); // check length ??

	double *stochNormModes = new double[N_Total_ModeDOF]();
	TFEVectFunct2D *VelocityStochNorm;
	VelocityStochNorm = new TFEVectFunct2D(VelocityMode_FeSpace,(char *)"U_Mode_StochNorm",(char *)"Stochastically normalized modes",stochNormModes,N_M,2*subDim);

	u1Mean = Velocity_Mean->GetComponent(0);
	u2Mean = Velocity_Mean->GetComponent(1);

	TFEVectFunct2D **VelocityModeAll = new TFEVectFunct2D *[subDim];
	for (int subD = 0; subD < subDim; subD++)
		VelocityModeAll[subD] = new TFEVectFunct2D(VelocityMode_FeSpace, (char *)"U_Mode", (char *)"Mode Component", solModeAll + (subD * 2 * N_M), N_M, 2);

	//======================================================================
	// /DO - SystemMatrix construction and solution
	//======================================================================
	// Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
	// Solver: AMG_SOLVE (or) GMG  (or) DIRECT
	TSystemTBE2D *SystemMatrix, *SystemMatrix_Mean, *SystemMatrix_Mode;
	SystemMatrix_Mean = new TSystemTBE2D(VelocityMean_FeSpace, Velocity_Mean, solMean, rhsMean, GALERKIN, DIRECT);

	TSystemTBE2D **SystemMatrixModeAll = new TSystemTBE2D *[subDim];
	for (int subD = 0; subD < subDim; subD++)
		SystemMatrixModeAll[subD] = new TSystemTBE2D(VelocityMode_FeSpace, VelocityModeAll[subD], solModeAll + (subD * 2 * N_M), rhsModeAll + (subD * 2 * N_M), GALERKIN, DIRECT);

	fefctMean[0] = u1Mean;
	fefctMean[1] = u2Mean;
	fespMean[0] = VelocityMean_FeSpace;

	BEaux_mean = new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2, TimeNSN_ParamFct2,
								 TimeNSN_FEValues2, fespMean, fefctMean, TimeNSFct2, TimeNSFEFctIndex2,
								 TimeNSFEMultiIndex2, TimeNSN_Params2, TimeNSBeginParam2);

	// aux for calculating the error
	// if (TDatabase::ParamDB->MEASURE_ERRORS)
	// {
	// 	BEaux_error_mean = new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
	// 									   TimeNSN_ParamFct2,
	// 									   TimeNSN_FEValues2,
	// 									   fespMean, FeFct,
	// 									   TimeNSFct2,
	// 									   TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
	// 									   TimeNSN_Params2, TimeNSBeginParam2);
	// }

	BEaux_modeAll = new TAuxParam2D *[subDim];
	fespMode[0] = VelocityMode_FeSpace;

	for (int subD = 0; subD < subDim; subD++)
	{
		fefctModeAll[subD][0] = VelocityModeAll[subD]->GetComponent(0); // changed ger from subDim to subD
		fefctModeAll[subD][1] = VelocityModeAll[subD]->GetComponent(1);

		BEaux_modeAll[subD] = new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2, TimeNSN_ParamFct2,
											  TimeNSN_FEValues2, fespMode, fefctModeAll[subD], TimeNSFct2, TimeNSFEFctIndex2,
											  TimeNSFEMultiIndex2, TimeNSN_Params2, TimeNSBeginParam2);
	}

	// aux for calculating the error
	// if (TDatabase::ParamDB->MEASURE_ERRORS)
	// {
	// 	BEaux_error_mode = new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
	// 									   TimeNSN_ParamFct2,
	// 									   TimeNSN_FEValues2,
	// 									   fespMode, FeFct,
	// 									   TimeNSFct2,
	// 									   TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
	// 									   TimeNSN_Params2, TimeNSBeginParam2);
	// }

	// SystemMatrix_Mode->Init(DO_Mode_Equation_Coefficients, BoundCondition, U1BoundValue, U2BoundValue, BEaux_mode, BEaux_error_mode);

	// -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
	//------------------------------------------ MEAN EQUATION SETUP -----------------------------------------------------//
	// initilize the system matrix with the functions defined in Example file
	SystemMatrix_Mean->Init(DO_Mean_Equation_Coefficients, BoundCondition, U1BoundValue, U2BoundValue, BEaux_mean, BEaux_error_mean);

	SystemMatrix_Mean->Assemble(solMean, rhsMean);

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
	for (int subD = 0; subD < subDim; subD++)
	{
		SystemMatrixModeAll[subD]->Init(DO_Mode_Equation_Coefficients, BoundCondition, U1BoundValue, U2BoundValue, BEaux_modeAll[subD], BEaux_error_mode);
		SystemMatrixModeAll[subD]->Assemble(solModeAll + (subD * 2 * N_M), rhsModeAll + (subD * 2 * N_M));
	}

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
	TDatabase::ParamDB->N_Subspace_Dim = subDim;
	TDatabase::ParamDB->COVARIANCE_MATRIX_DO = new double[subDim * subDim]();
	TDatabase::ParamDB->COSKEWNESS_MATRIX_DO = new double[subDim * subDim * subDim]();
	TDatabase::ParamDB->COVARIANCE_INVERSE_DO = new double[subDim * subDim]();

	CalcCovarianceMatx(CoeffVector);

	CalcCoskewnessMatx(CoeffVector);
	InvertCov();

	// 	// double* ModeVector_OldRHS = new double[N_U]();

	// // Set up a FE VECT FUNCTION TO STORE ALL THE Components of CTilde

	//     // TFEVectFunct2D* linModesFeVectFunct =

	// int TimeLinear_FESpaces_DO = 1;
	// int TimeLinear_Fct_DO = 1; // \tilde(C)
	// int TimeLinear_ParamFct_DO = 1;
	// int TimeLinear_FEValues_DO = 3;
	// int TimeLinear_Params_DO = 3;
	// int TimeNSFEFctIndex_DO[3] = {0, 0, 0};
	// MultiIndex2D TimeNSFEMultiIndex_DO[3] = {D00, D01, D10};
	// ParamFct *TimeNSFct_DO[1] = {DO_Mode_RHS_Aux_Param};
	// int TimeNSBeginParam_DO[1] = {0};

	// TFEFunction2D *fefct_RHS[4];
	// TFESpace2D *fesp_RHS[2];

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

	TOutput2D **OutputModeAll = new TOutput2D *[subDim];

	for (int sd = 0; sd < subDim; sd++)
	{
		OutputModeAll[sd] = new TOutput2D(2, 2, 1, 1, Domain);
		OutputModeAll[sd]->AddFEVectFunct(VelocityModeAll[sd]);
	}

	// OutputMode = new TOutput2D(2, 2*subDim, 1, 1, Domain);
	// OutputMode->AddFEVectFunct(Velocity_Mode);

	int meanimg = 0;
	int *modeimg = new int[subDim]();
	int coeffimg = 0;

	m = 0;

	std::string meanBaseName = "Mean/Mean_NRealisations_";
	std::string modeBaseName = "Modes/Mode_NRealisations_";
	std::string coeffBaseName = "Coefficients/Coeff_NRealisations_";
	std::string IPMeanBaseName = "IPMatrices/IPMean_NRealisations_";
	std::string IPModeBaseName = "IPMatrices/IPMode_NRealisations_";

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
	fileoutMean = generateFileName(meanBaseName, m, N_Realisations);
	printToTxt(fileoutMean, solMean, 2 * N_U, 1, 'C');

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
	fileoutMode = generateFileName(modeBaseName, m, N_Realisations);
	printToTxt(fileoutMode, solModeAll, 2 * N_U, subDim, 'C');

	// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	// coll->GetHminHmax(&hmin, &hmax);
	// OutPut("h_min : " << hmin << " h_max : " << hmax << endl);

	//======================================================================
	// time disc loop
	//======================================================================
	// parameters for time stepping scheme

	fileoutCoeff = generateFileName(coeffBaseName, m, N_Realisations);
	printToTxt(fileoutCoeff, CoeffVector, N_Realisations, subDim, 'C');

	double *u1_Mode = new double[N_U * subDim]();
	double *u2_Mode = new double[N_U * subDim]();
	for (int i = 0; i < N_U; i++)
	{
		for (int j = 0; j < subDim; j++)
		{
			u1_Mode[j * N_U + i] = solModeAll[(2 * j * N_U) + i]; //??
			u2_Mode[j * N_U + i] = 0;							  //??
		}
	}

	double *u1_Mean = new double[N_U * 1]();
	double *u2_Mean = new double[N_U * 1]();

	for (int i = 0; i < N_U; i++)
		u1_Mean[i] = solMean[i];

	std::string fileoutIPMode;
	std::string fileoutIPMean;

	fileoutIPMean = generateFileName(IPMeanBaseName, m, N_Realisations);
	fileoutIPMode = generateFileName(IPModeBaseName, m, N_Realisations);

	calcIPMatx(IPMatxMode, u1_Mode, N_U, subDim, 'C');
	calcIPMatx(IPMatxMean, u1_Mean, N_U, 1, 'C');

	printToTxt(fileoutIPMean, IPMatxMean, N_M, 1, 'C');
	printToTxt(fileoutIPMode, IPMatxMode, N_M, subDim, 'C');

	N_SubSteps = GetN_SubSteps();
	oldtau = 1.;
	end_time = TDatabase::TimeDB->ENDTIME;
	limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
	Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
	memset(AllErrors, 0, 7. * SizeOfDouble);

	while (TDatabase::TimeDB->CURRENTTIME < end_time)
	{ // time cycle
		// fileoutCoeff = generateFileName(coeffBaseName, m, N_Realisations);

		// printToTxt(fileoutCoeff, CoeffVector, N_Realisations, subDim, 'C');

		// fileoutMode = generateFileName(modeBaseName, m, N_Realisations);

		// printToTxt(fileoutMode, solModeAll, 2 * N_U, subDim, 'C');

		// fileoutMean = generateFileName(meanBaseName, m, N_Realisations);

		// printToTxt(fileoutMean, solMean, 2 * N_U, 1, 'C');
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

			OutPut(endl << "CURRENT TIME: ");
			OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

			// copy sol, rhs to olssol, oldrhs
			memcpy(old_rhsMean, rhsMean, N_Total_MeanDOF * SizeOfDouble);
			memcpy(old_solMean, solMean, N_Total_MeanDOF * SizeOfDouble);
			normalizeStochasticModes(VelocityMode_FeSpace,Velocity_Mode,subDim,stochNormModes);
			DO_Mean_RHS(VelocityMean_FeSpace, Velocity_Mode, subDim, rhsMean, N_U);

			// assemble only rhs, nonlinear matrix for NSE will be assemble in fixed point iteration
			SystemMatrix_Mean->Assemble(solMean, rhsMean);

			// scale B matices and assemble NSE-rhs based on the \theta time stepping scheme
			//  SystemMatrix_Mean->AssembleSystMat(oldrhs, rhs, sol);
			SystemMatrix_Mean->AssembleSystMat(old_rhsMean, rhsMean, solMean);
			oldtau = tau;

			// calculate the residual
			memset(defect, 0, N_Total_MeanDOF * SizeOfDouble);
			// SystemMatrix_Mean->GetTBEResidual(sol, defect);
			SystemMatrix_Mean->GetTBEResidual(solMean, defect);

			residual = Ddot(N_Total_MeanDOF, defect, defect);
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
				memset(defect, 0, N_Total_MeanDOF * SizeOfDouble);
				//  SystemMatrix_Mean->GetTBEResidual(sol, defect);
				SystemMatrix_Mean->GetTBEResidual(solMean, defect);

				residual = Ddot(N_Total_MeanDOF, defect, defect);
				OutPut("nonlinear iteration step " << setw(3) << j);
				OutPut(setw(14) << sqrt(residual) << endl);

				if (sqrt(residual) <= limit)
					break;

			} // for(j=1;j<=Max_It;j++)

			// restore the mass matrix for the next time step
			// SystemMatrix_Mean->RestoreMassMat();
			SystemMatrix_Mean->RestoreMassMat();

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

		// OutPut("Ddot Coeff Vector Before Coeff Solve:" << Ddot(N_Realisations*subDim, CoeffVector, CoeffVector) << endl);
		for (int subSpaceNum = 0; subSpaceNum < subDim; subSpaceNum++)
		{

			DO_CoEfficient(Velocity_FeSpace, Velocity_Mode, FeVector_Coefficient, Velocity_Mean, subDim, subSpaceNum, N_Realisations);
		}

		// OutPut("Ddot Coeff Vector After Coeff Solve:" << Ddot(N_Realisations*subDim, CoeffVector, CoeffVector) << endl);
		CalcCovarianceMatx(CoeffVector);
		CalcCoskewnessMatx(CoeffVector);
		InvertCov();
		for (int subSpaceNum = 0; subSpaceNum < subDim; subSpaceNum++)
		{
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
				// TDatabase::TimeDB->CURRENTTIME += tau;

				// subspace loop
				double *modeSolution_i = solModeAll + subSpaceNum * 2 * N_M; // This works for column major
				double *modeSolution_rhs = rhsModeAll + subSpaceNum * 2 * N_M;

				// copy sol, rhs to olssol, oldrhs
				memcpy(old_rhsMode, modeSolution_rhs, 2 * N_M * SizeOfDouble);
				// memcpy(old_solMode, modeSolution_i, N_Total_MeanDOF * SizeOfDouble);

				// Assemble rhs
				DO_Mode_RHS(VelocityMode_FeSpace, Velocity_Mean, Velocity_Mode, subDim, modeSolution_rhs, subSpaceNum);

				SystemMatrixModeAll[subSpaceNum]->Assemble(modeSolution_i, modeSolution_rhs);
				//   }

				// scale B matices and assemble NSE-rhs based on the \theta time stepping scheme
				//  SystemMatrix_Mean->AssembleSystMat(oldrhs, rhs, sol);
				SystemMatrixModeAll[subSpaceNum]->AssembleSystMat(old_rhsMode, modeSolution_rhs, modeSolution_i);
				oldtau = tau;

				// calculate the residual
				memset(defect, 0, N_Total_MeanDOF * SizeOfDouble);
				// SystemMatrix_Mean->GetTBEResidual(sol, defect);
				SystemMatrixModeAll[subSpaceNum]->GetTBEResidual(modeSolution_i, defect);

				residual = Ddot(N_Total_MeanDOF, defect, defect);
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
					SystemMatrixModeAll[subSpaceNum]->Solve(modeSolution_i);

					// restore the mass matrix for the next nonlinear iteration
					//  SystemMatrix_Mean->RestoreMassMat();
					SystemMatrixModeAll[subSpaceNum]->RestoreMassMat();

					// assemble the system matrix with given aux, sol and rhs
					//  SystemMatrix_Mean->AssembleANonLinear(sol, rhs);
					SystemMatrixModeAll[subSpaceNum]->AssembleANonLinear(modeSolution_i, modeSolution_rhs);

					// assemble system mat, S = M + dt\theta_1*A
					//  SystemMatrix_Mean->AssembleSystMatNonLinear();
					SystemMatrixModeAll[subSpaceNum]->AssembleSystMatNonLinear();

					// get the residual
					memset(defect, 0, N_Total_MeanDOF * SizeOfDouble);
					//  SystemMatrix_Mean->GetTBEResidual(sol, defect);
					SystemMatrixModeAll[subSpaceNum]->GetTBEResidual(modeSolution_i, defect);

					residual = Ddot(N_Total_MeanDOF, defect, defect);
					//  OutPut("nonlinear iteration step " << setw(3) << j);
					//  OutPut(setw(14) << sqrt(residual) << endl);

					if (sqrt(residual) <= limit)
						break;

				} // for(j=1;j<=Max_It;j++)
				SystemMatrixModeAll[subSpaceNum]->RestoreMassMat();
			} // l loop end

		} //  subspace loop end
		//======================================================================
		// produce outout
		//======================================================================
		if (m == 1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
		{
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
		}

		// xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
		// reorthonormalizeC(solModeAll,N_Total_ModeDOF,subDim);
        memcpy(solMode, solModeAll, N_Total_ModeDOF * SizeOfDouble);

	} // while(TDatabase::TimeDB->CURRENTTIME< e

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

	fileoutCoeff = generateFileName(coeffBaseName, m, N_Realisations);

	printToTxt(fileoutCoeff, CoeffVector, N_Realisations, subDim, 'C');

	fileoutMode = generateFileName(modeBaseName, m, N_Realisations);

	printToTxt(fileoutMode, solModeAll, 2 * N_U, subDim, 'C');

	fileoutMean = generateFileName(meanBaseName, m, N_Realisations);

	printToTxt(fileoutMean, solMean, 2 * N_U, 1, 'C');
	TDatabase::TimeDB->CURRENTTIME = 0;

	CloseFiles();
	std::string PyInFile = "PyIn.txt";

	std::ofstream PyIn;
	PyIn.open(PyInFile);
	PyIn << N_Realisations << endl
		 << subDim << endl
		 << m << endl;

	double *RealizationVectorCopy = new double[N_U * N_Realisations]();
	memcpy(RealizationVectorCopy, RealizationVector, N_U * N_Realisations * SizeOfDouble);

	double *MeanVectorMC = new double[N_U * 1]();
	calcMeanRealization(RealizationVectorCopy, MeanVectorMC, N_Realisations, N_U);

	double *stdDevVectorMC = new double[N_U * 1]();
	calcStdDevRealization(RealizationVectorCopy, stdDevVectorMC, N_Realisations, N_U);

	double *CompositeVectorMC = new double[N_U * 7]();
	for (int i = 0; i < N_U; i++)
	{
		CompositeVectorMC[i] = MeanVectorMC[i];
		CompositeVectorMC[N_U + i] = MeanVectorMC[i] + stdDevVectorMC[i];
		CompositeVectorMC[2 * N_U + i] = MeanVectorMC[i] - stdDevVectorMC[i];
		CompositeVectorMC[3 * N_U + i] = MeanVectorMC[i] + (2 * stdDevVectorMC[i]);
		CompositeVectorMC[4 * N_U + i] = MeanVectorMC[i] - (2 * stdDevVectorMC[i]);
		CompositeVectorMC[5 * N_U + i] = MeanVectorMC[i] + (3 * stdDevVectorMC[i]);
		CompositeVectorMC[6 * N_U + i] = MeanVectorMC[i] - (3 * stdDevVectorMC[i]);
	}

	double *solMC = new double[2 * N_U]();
	double *rhsMC = new double[2 * N_U]();
	double *oldrhsMC = new double[2 * N_U]();
	double *defectMC = new double[2 * N_U]();

	TFEVectFunct2D *VelocityMC = new TFEVectFunct2D(Velocity_FeSpace, (char *)"MCSol", (char *)"Solution of Monte Carlo 7", solMC, N_U, 2);
	TFEFunction2D *u1MC = VelocityMC->GetComponent(0);
	TFEFunction2D *u2MC = VelocityMC->GetComponent(1);

	TSystemTBE2D *SystemMatrixMC = new TSystemTBE2D(Velocity_FeSpace, VelocityMC, solMC, rhsMC, Disctype, DIRECT);

	FeFct[0] = u1MC;
	FeFct[1] = u2MC;
	fesp[0] = Velocity_FeSpace;
	BEaux = new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2, TimeNSN_ParamFct2,
							TimeNSN_FEValues2, fesp, FeFct, TimeNSFct2, TimeNSFEFctIndex2,
							TimeNSFEMultiIndex2, TimeNSN_Params2, TimeNSBeginParam2);

	// aux for calculating the error
	if (TDatabase::ParamDB->MEASURE_ERRORS)
	{
		BEaux_error = new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
									  TimeNSN_ParamFct2,
									  TimeNSN_FEValues2,
									  fesp, FeFct,
									  TimeNSFct2,
									  TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
									  TimeNSN_Params2, TimeNSBeginParam2);
	}

	SystemMatrixMC->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue, BEaux, BEaux_error);

	TOutput2D *OutputMC = new TOutput2D(2, 2, 1, 1, Domain);
	OutputMC->AddFEVectFunct(VelocityMC);

	int N_Composite = 7;
	int *imgm = new int[N_Composite]();
	for (int RealNo = 0; RealNo < 7; RealNo++)
	{
		std::string filename;
		cout << " ============================================================================================================= " << endl;
		switch (RealNo)
		{
		case 0:
			cout << "Solving for Mean Solution" << endl;
			filename = "MonteCarlo_Mean_NR_";
			break;
		case 1:
			cout << "Solving for Mean + sigma Solution" << endl;
			filename = "MonteCarlo_MeanPlusSigma_NR_";
			break;
		case 2:
			cout << "Solving for Mean - sigma Solution" << endl;
			filename = "MonteCarlo_MeanMinusSigma_NR_";
			break;
		case 3:
			cout << "Solving for Mean + 2*sigma Solution" << endl;
			filename = "MonteCarlo_MeanPlus2Sigma_NR_";
			break;
		case 4:
			cout << "Solving for Mean - 2*sigma Solution" << endl;
			filename = "MonteCarlo_MeanMinus2Sigma_NR_";
			break;
		case 5:
			cout << "Solving for Mean + 3*sigma Solution" << endl;
			filename = "MonteCarlo_MeanPlus3Sigma_NR_";
			break;
		case 6:
			cout << "Solving for Mean - 3*sigma Solution" << endl;
			filename = "MonteCarlo_MeanMinus3Sigma_NR_";
			break;
		}
		cout << " ============================================================================================================= " << endl;

		VtkBaseName = const_cast<char *>(filename.c_str());

		for (int i = 0; i < N_U; i++)
			solMC[i] = CompositeVectorMC[i + RealNo * N_U];

		os.seekp(std::ios::beg);
		if (imgm[RealNo] < 10)
			os << "VTK"
			   << "/" << VtkBaseName << ".0000" << imgm[RealNo] << ".vtk" << ends;
		else if (imgm[RealNo] < 100)
			os << "VTK"
			   << "/" << VtkBaseName << ".000" << imgm[RealNo] << ".vtk" << ends;
		else if (imgm[RealNo] < 1000)
			os << "VTK"
			   << "/" << VtkBaseName << ".00" << imgm[RealNo] << ".vtk" << ends;
		else if (imgm[RealNo] < 10000)
			os << "VTK"
			   << "/" << VtkBaseName << ".0" << imgm[RealNo] << ".vtk" << ends;
		else
			os << "VTK"
			   << "/" << VtkBaseName << "." << imgm[RealNo] << ".vtk" << ends;
		OutputMC->WriteVtk(os.str().c_str());

		fileoutMC = generateFileName("MonteCarlo/" + filename, imgm[RealNo], N_Realisations);
		printToTxt(fileoutMC, solMC, N_U, 1, 'C');
		imgm[RealNo]++;

		SystemMatrixMC->Assemble(solMC, rhsMC);

		m = 0;
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
				memcpy(oldrhsMC, rhsMC, 2 * N_U * SizeOfDouble);

				// assemble only rhs, nonlinear matrix for NSE will be assemble in fixed point iteration
				// not needed if rhs is not time-dependent
				if (m != 1)
				{
					//  SystemMatrix_Mean->AssembleA();
					SystemMatrixMC->AssembleA();
				}
				else
				{
					//  SystemMatrix_Mean->Assemble(sol, rhs);
					SystemMatrixMC->Assemble(solMC, rhsMC);
				}

				// scale B matices and assemble NSE-rhsMC based on the \theta time stepping scheme
				//  SystemMatrix_Mean->AssembleSystMat(oldrhs, rhsMC, solMC);
				SystemMatrixMC->AssembleSystMat(oldrhsMC, rhsMC, solMC);
				oldtau = tau;

				// calculate the residual
				// SystemMatrix_Mean->GetTBEResidual(solMC, defect);
				SystemMatrixMC->GetTBEResidual(solMC, defectMC);

				residual = Ddot(2 * N_U, defectMC, defectMC);
				// OutPut("Nonlinear iteration step   0");
				// OutPut(setw(14) << sqrt(residual) << endl);

				//======================================================================
				// Solve the system
				// Nonlinear iteration of fixed point type
				//======================================================================
				for (j = 1; j <= Max_It; j++)
				{
					// Solve the NSE system
					//  SystemMatrix_Mean->Solve(solMC);
					SystemMatrixMC->Solve(solMC);

					// restore the mass matrix for the next nonlinear iteration
					//  SystemMatrix_Mean->RestoreMassMat();
					SystemMatrixMC->RestoreMassMat();

					// assemble the system matrix with given aux, solMC and rhsMC
					//  SystemMatrix_Mean->AssembleANonLinear(solMC, rhsMC);
					SystemMatrixMC->AssembleANonLinear(solMC, rhsMC);

					// assemble system mat, S = M + dt\theta_1*A
					//  SystemMatrix_Mean->AssembleSystMatNonLinear();
					SystemMatrixMC->AssembleSystMatNonLinear();

					// get the residual
					memset(defectMC, 0, 2 * N_U * SizeOfDouble);
					//  SystemMatrix_Mean->GetTBEResidual(solMC, defect);
					SystemMatrixMC->GetTBEResidual(solMC, defectMC);

					residual = Ddot(2 * N_U, defectMC, defectMC);
					//  OutPut("nonlinear iteration step " << setw(3) << j);
					//  OutPut(setw(14) << sqrt(residual) << endl);

					if (sqrt(residual) <= limit)
						break;

				} // for(j=1;j<=Max_It;j++)

				// restore the mass matrix for the next time step
				// SystemMatrix_Mean->RestoreMassMat();
				SystemMatrixMC->RestoreMassMat();

			} // for(l=0;l<N_SubSteps;
			  //======================================================================
			  // measure errors to known solution
			  //======================================================================
			if (TDatabase::ParamDB->MEASURE_ERRORS)
			{
				// SystemMatrix_Mean->MeasureErrors(ExactU1, ExactU2, AllErrors);
				SystemMatrixMC->MeasureErrors(ExactU1, ExactU2, AllErrors);

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
					if (imgm[RealNo] < 10)
						os << "VTK"
						   << "/" << VtkBaseName << ".0000" << imgm[RealNo] << ".vtk" << ends;
					else if (imgm[RealNo] < 100)
						os << "VTK"
						   << "/" << VtkBaseName << ".000" << imgm[RealNo] << ".vtk" << ends;
					else if (imgm[RealNo] < 1000)
						os << "VTK"
						   << "/" << VtkBaseName << ".00" << imgm[RealNo] << ".vtk" << ends;
					else if (imgm[RealNo] < 10000)
						os << "VTK"
						   << "/" << VtkBaseName << ".0" << imgm[RealNo] << ".vtk" << ends;
					else
						os << "VTK"
						   << "/" << VtkBaseName << "." << imgm[RealNo] << ".vtk" << ends;
					OutputMC->WriteVtk(os.str().c_str());

					fileoutMC = generateFileName("MonteCarlo/" + filename, imgm[RealNo], N_Realisations);
					printToTxt(fileoutMC, solMC, N_U, 1, 'C');
					imgm[RealNo]++;
				}

		} // while(TDatabase::TimeDB->CURRENTTIME< e
		if (TDatabase::ParamDB->WRITE_VTK)
		{
			os.seekp(std::ios::beg);
			if (imgm[RealNo] < 10)
				os << "VTK"
				   << "/" << VtkBaseName << ".0000" << imgm[RealNo] << ".vtk" << ends;
			else if (imgm[RealNo] < 100)
				os << "VTK"
				   << "/" << VtkBaseName << ".000" << imgm[RealNo] << ".vtk" << ends;
			else if (imgm[RealNo] < 1000)
				os << "VTK"
				   << "/" << VtkBaseName << ".00" << imgm[RealNo] << ".vtk" << ends;
			else if (imgm[RealNo] < 10000)
				os << "VTK"
				   << "/" << VtkBaseName << ".0" << imgm[RealNo] << ".vtk" << ends;
			else
				os << "VTK"
				   << "/" << VtkBaseName << "." << imgm[RealNo] << ".vtk" << ends;
			OutputMC->WriteVtk(os.str().c_str());

			fileoutMC = generateFileName("MonteCarlo/" + filename, imgm[RealNo], N_Realisations);
			printToTxt(fileoutMC, solMC, N_U, 1, 'C');
			imgm[RealNo]++;
		}
	} // Realization Loop
	TDatabase::TimeDB->CURRENTTIME = 0;

	delete[] RealizationVectorCopy;
	delete[] MeanVectorMC;
	delete[] stdDevVectorMC;
	delete[] CompositeVectorMC;

	CloseFiles();

	return 0;

} // end main
