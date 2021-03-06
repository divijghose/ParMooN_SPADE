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

	double *sol, *solold, *rhs, *oldrhs, *defect, t1, t2, residual, impuls_residual;
	double limit, AllErrors[7], end_time, oldtau, tau;

	TDomain *Domain;
	TDatabase *Database = new TDatabase();
	TFEDatabase2D *FEDatabase = new TFEDatabase2D();
	TCollection *coll, *mortarcoll = NULL;
	TFESpace2D *Velocity_FeSpace, *Pressure_FeSpace, *fesp[2];
	TFEVectFunct2D *Velocity;
	TFEFunction2D *u1, *u2, *Pressure, *fefct[2];
	TOutput2D *Output;
	TOutput2D *meanOutput;
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
    Disctype = 1;
	coll = Domain->GetCollection(It_Finest, 0);
	N_Cells = coll->GetN_Cells();
	OutPut("N_Cells : " << N_Cells << endl);
	std::ofstream pyInputFile;
	pyInputFile.open("pyReadIn.txt");

	// fespaces for velocity and pressure
	GetVelocityAndPressureSpace(coll, BoundCondition, mortarcoll, Velocity_FeSpace,
								Pressure_FeSpace, &pressure_space_code,
								TDatabase::ParamDB->VELOCITY_SPACE,
								TDatabase::ParamDB->PRESSURE_SPACE);

	// defaulty inf-sup pressure space will be selected based on the velocity space, so update it in database
	TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
	velocity_space_code = TDatabase::ParamDB->VELOCITY_SPACE;

	N_U = Velocity_FeSpace->GetN_DegreesOfFreedom();
	pyInputFile << N_U << endl;
	N_P = Pressure_FeSpace->GetN_DegreesOfFreedom();
	pyInputFile << N_P << endl;

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
	solold = new double[N_TotalDOF];
	rhs = new double[N_TotalDOF];
	oldrhs = new double[N_TotalDOF];

	memset(sol, 0, N_TotalDOF * SizeOfDouble);
	memset(solold, 100, N_TotalDOF * SizeOfDouble);
	memset(rhs, 0, N_TotalDOF * SizeOfDouble);

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

	// Looping Starts here 
	int N_Realisations = 5;
	pyInputFile << N_Realisations << endl;
	std::ofstream fileoutSolution;
	fileoutSolution.open("FinalSolution_AllRealisations.txt");
    srand(time(0));
    std::ofstream fileout;
    std::ofstream reNR;
	reNR.open("ReNo.txt");

	std::ofstream fileoutTimeStep;

	int* RENR = new int[N_Realisations];
	for ( int i = 0 ; i < N_Realisations; i++)
		RENR[i] = rand() % 2000;
	
	for ( int Realisation = 0 ; Realisation < N_Realisations ; Realisation++)
	{   

        TDatabase::ParamDB->RE_NR =   RENR[Realisation];
        cout << " ------------------------------------- REALISATION ----------------------------- " <<endl;
        cout << " ----------------------" <<  " NO : " << Realisation << "  RE_NR : " <<  TDatabase::ParamDB->RE_NR << " ---------------------- " <<endl;
		// Ste variables as Zero
		memset(sol,0,SizeOfDouble*N_TotalDOF);
		memset(solold,100,SizeOfDouble*N_TotalDOF);
		memset(rhs,0,SizeOfDouble*N_TotalDOF);

		std::string name = "ReNr_" + std::to_string(int(TDatabase::ParamDB->RE_NR)) ;
		mkdir(name.c_str(), 0777);

        std::string a = std::to_string(Realisation);
        reNR << a << " " << TDatabase::ParamDB->RE_NR << endl;
		// Change Reynolds number for the iteration

		SystemMatrix->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue, aux, NSEaux_error);

		SystemMatrix->Assemble(sol, rhs);


		// assemble M, A matrices and rhs
		SystemMatrix->Assemble(sol, rhs);
		

		//======================================================================
		// produce outout
		//======================================================================
		VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
		std::string filename = "ReNr_" + std::to_string(int(TDatabase::ParamDB->RE_NR)) ;
		VtkBaseName   = const_cast<char*>(filename.c_str());
		Output = new TOutput2D(2, 2, 1, 1, Domain);
		meanOutput = new TOutput2D(2, 2, 1, 1, Domain);
		img = 0;


		
		Output->AddFEVectFunct(Velocity);
		meanOutput->AddFEVectFunct(Velocity);
		Output->AddFEFunction(Pressure);
		meanOutput->AddFEFunction(Pressure);
		const char* outputfilename;
		const char* meanoutputfilename;

		if (TDatabase::ParamDB->WRITE_VTK)
		{
			os.seekp(std::ios::beg);
			if (img < 10)
				os << name.c_str()  << "/"  << VtkBaseName << ".0000" << img << ".vtk" << ends;
			else if (img < 100)
				os << "VTK/" << VtkBaseName << ".000" << img << ".vtk" << ends;
			else if (img < 1000)
				os << name.c_str()  << "/" << VtkBaseName << ".00" << img << ".vtk" << ends;
			else if (img < 10000)
				os << name.c_str()  << "/" << VtkBaseName << ".0" << img << ".vtk" << ends;
			else
				os << name.c_str()  << "/" << VtkBaseName << "." << img << ".vtk" << ends;
			outputfilename = os.str().c_str();
			Output->WriteVtk(os.str().c_str());
			
			img++;
		}

		
		fileout.open( a + ".txt");
		for( int i=0 ; i < N_TotalDOF ; i++)  fileout<<sol[i] << "\t";
		fileout<<endl;

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
		double norm = 100;
		// time loop starts
		while (TDatabase::TimeDB->CURRENTTIME < end_time)
		{ // time cycle

			
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

				// OutPut(endl
				// 	<< "CURRENT TIME: ");
				// OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

				//copy sol, rhs to olssol, oldrhs
				memcpy(oldrhs, rhs, N_TotalDOF * SizeOfDouble);

				// assemble only rhs, nonlinear matrix for NSE will be assemble in fixed point iteration
				// not needed if rhs is not time-dependent
				// if (m != 1)
				// {
					SystemMatrix->AssembleRhs(sol, rhs);
				// }
				// else
				// {
					SystemMatrix->Assemble(sol, rhs);
				// }

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
				// OutPut("Nonlinear iteration step   0");
				// OutPut(setw(14) << impuls_residual);
				// OutPut(setw(14) << residual - impuls_residual);
				// OutPut(setw(14) << sqrt(residual) << endl);

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
				
				// cout << "Solution Norm = " << Ddot(N_TotalDOF, sol, sol) << endl;

			} // for(l=0;l<N_SubSteps;
			//======================================================================
			// measure errors to known solution
			//======================================================================


			//======================================================================
			// produce outout
			//======================================================================
			if (m == 1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
				if (TDatabase::ParamDB->WRITE_VTK)
				{
					os.seekp(std::ios::beg);
					if (img < 10)
						os << name.c_str()  << "/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
					else if (img < 100)
						os << name.c_str()  << "/" << VtkBaseName << ".000" << img << ".vtk" << ends;
					else if (img < 1000)
						os << name.c_str()  << "/" << VtkBaseName << ".00" << img << ".vtk" << ends;
					else if (img < 10000)
						os << name.c_str()  << "/" << VtkBaseName << ".0" << img << ".vtk" << ends;
					else
						os << name.c_str()  << "/" << VtkBaseName << "." << img << ".vtk" << ends;
					outputfilename = os.str().c_str();
					Output->WriteVtk(os.str().c_str());
					img++;
				}

				for( int i=0 ; i < N_TotalDOF ; i++)  fileout<<sol[i] << "\t";
				fileout<<endl;



		} // while(TDatabase::TimeDB->CURRENTTIME< e

		//======================================================================
		// produce final outout
		//======================================================================
		if (TDatabase::ParamDB->WRITE_VTK)
		{
			os.seekp(std::ios::beg);
			if (img < 10)
				os << name.c_str()  << "/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
			else if (img < 100)
				os << name.c_str()  << "/" << VtkBaseName << ".000" << img << ".vtk" << ends;
			else if (img < 1000)
				os << name.c_str()  << "/" << VtkBaseName << ".00" << img << ".vtk" << ends;
			else if (img < 10000)
				os << name.c_str()  << "/" << VtkBaseName << ".0" << img << ".vtk" << ends;
			else
				os << name.c_str()  << "/" << VtkBaseName << "." << img << ".vtk" << ends;
			Output->WriteVtk(os.str().c_str());
			img++;
		}

		double solnorm = 0.0;
		solnorm = Ddot(N_TotalDOF, sol, sol);
		cout << "Solution Norm = " << solnorm << endl;

        for( int i=0 ; i < N_TotalDOF ; i++)  fileout<<sol[i] << "\t";
        fileout<<endl;
        fileout.close();

		for ( int i=0 ; i < N_TotalDOF;i++)  fileoutSolution << sol[i] << ",";
		fileoutSolution<<endl;

		


        TDatabase::TimeDB->CURRENTTIME = 0;

	}  // END OF REALISATION LOOP 
	pyInputFile << m << endl;
	pyInputFile.close();
	fileoutSolution.close();
	
    reNR.close();
   

	CloseFiles();

	return 0;
} // end main
