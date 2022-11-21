/**
 * @file TBE2D_ParMooN_MC.C
 * @brief Purpose:  Main program for solving Monte Carlo runs of time-dependent Burger's equation in ParMooN.			   
 * @authors Sashikumaar Ganesan
 * @authors Divij Ghose
 * @authors Thivin Anandh
 * @bug No known bugs
 */

// ===========================================================================//
//  Monte Carlo Solution of Burgers' Equation //
// ===========================================================================//

// =======================================================================
//
// Purpose:     Main program for solving Monte Carlo runs of time-dependent Burger's equation in ParMooN.	
//
// Authors:      Sashikumaar Ganesan, Thivin Anandh, Divij Ghose
//
// History:     1> First iteration implemented on 01.04.2022
//              2> Sine Wave standard deviation function implemented on 12.04.2022
//              
//		        
//
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
#include "../Examples/Monte_Carlo/burgers_mc.h" // smooth sol in unit square

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  // ======================================================================
  //  declaration of variables
  // ======================================================================
  int i, j, l, m, N_Cells, ORDER, N_U, N_M, N_L, N_TotalDOF, img=1, N_Modes, N_Realiz;
  int N_Total_MeanDOF, Max_It, NSEType, N_SubSteps, Disctype;
  
  double *sol, *sol_mode, *rhs, *oldrhs, *defect, t1, t2, residual, impuls_residual;
  double limit, AllErrors[7], end_time, oldtau, tau;
  
  TDomain *Domain;
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D(); 
  TCollection *coll, *mortarcoll = NULL;
  // TFESpace2D *Velocity_FeSpace, *VelocityMode_FeSpace, *fesp[2], *fesp_mode[2];
  TFESpace2D *Velocity_FeSpace, *fesp[2];
  // TFEVectFunct2D *Velocity_Mean, *Velocity_Mode;
  TFEVectFunct2D *Velocity;
  // TFEFunction2D *u1_mean, *u2_mean, *fefct[2];
  TFEFunction2D *u1, *u2, *fefct[2];
  TOutput2D *Output;
  // TSystemTBE2D *SystemMatrix_Mean;
  TSystemTBE2D *SystemMatrix;
  // TSystemTBE_Mode2D *SystemMatrix_Mode;
  // TFEFunction2D *FeFct[2], *FeFct_Mode[4];
  TFEFunction2D *FeFct[2];
  TAuxParam2D *BEaux, *BE_Modeaux, *BEaux_error;

  const char vtkdir[] = "VTK"; 
  char *PsBaseName, *VtkBaseName, *GEO;
  // char UString[] = "u_mean";
  char UString[] = "usol";
  char UMString[] = "u_mode";
  char NameString[] = "UQ";
  
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
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();  
  
  if(TDatabase::ParamDB->WRITE_PS)
   {
    // write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
   }
   
//=========================================================================
// construct all finite element spaces
//=========================================================================
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  Disctype = TDatabase::ParamDB->DISCTYPE;

  // N_Modes = TDatabase::ParamDB->P10;
  // N_Realiz = TDatabase::ParamDB->P11;

  // VelocityMode = new TFEVectFunct2D*[N_Modes];

  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  OutPut("N_Cells : " << N_Cells <<endl);
  
  Velocity_FeSpace = new TFESpace2D(coll, (char*)"sol", (char*)"sol", BoundCondition, ORDER, NULL);
     
  // VelocityMode_FeSpace = new TFESpace2D(coll, (char*)"Mode", (char*)"Mode", BoundCondition, 1, NULL);

  

  N_U = Velocity_FeSpace->GetN_DegreesOfFreedom(); 
  // N_M = VelocityMode_FeSpace->GetN_DegreesOfFreedom(); 
  // N_Total_MeanDOF = 2*N_U;
  // N_TotalDOF = 2*(N_U + N_Modes*N_M);
  N_TotalDOF = 2*N_U;
  // OutPut("Dof Mean Velocity : "<< setw(10) << 2* N_U << endl);
  // OutPut("Dof Mode Velocity : "<< setw(10) << 2* N_M << endl);
  // OutPut("Total Dof all: "<< setw(10) << N_TotalDOF  << endl);
   OutPut("DOF : " << setw(10) << 2*N_U<<endl); 
//======================================================================
// construct all finite element functions
//======================================================================
    sol = new double[N_TotalDOF]; 
    rhs = new double[N_TotalDOF];
    oldrhs = new double[N_TotalDOF];
    defect = new double[N_TotalDOF];

    memset(sol, 0, N_TotalDOF*SizeOfDouble);
    memset(rhs, 0, N_TotalDOF*SizeOfDouble);

    Velocity = new TFEVectFunct2D(Velocity_FeSpace, UString,  UString,  sol, N_U, 2);
    u1 = Velocity->GetComponent(0);
    u2 = Velocity->GetComponent(1);

    u1->Interpolate(InitialU1);
    u2->Interpolate(InitialU2);

    /** mean velo */
    // Velocity_Mean = new TFEVectFunct2D(Velocity_FeSpace, UString,  UString,  sol, N_U, 2);
    // u1_mean = Velocity_Mean->GetComponent(0);
    // u2_mean = Velocity_Mean->GetComponent(1);

    //interpolate the initial solution
    // u1_mean->Interpolate(InitialU1Mean);
    // u2_mean->Interpolate(InitialU2Mean);

    /** mode velo */
    // sol_mode = sol+2*N_U;    
    // Velocity_Mode = new TFEVectFunct2D(VelocityMode_FeSpace, UString,  UString,  sol_mode, N_M, 2*N_Modes);

    // Sol_DO_Coeff = new double[N_Realiz*N_modes]


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
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
                double sig_r1 = exp ( - pow( ( 2*actual_x - 1 - disp),power)  / (E) )  / (2*3.14159265359 * sqrt(E))  * exp ( - pow(( 2*actual_y - 1-disp),power)  / (E) )  / (2*3.14159265359 * sqrt(E)) ;
                double sig_r2 = exp ( - pow(( 2*local_x - 1 -disp),power)  / (E) )  / (2*3.14159265359 * sqrt(E))  * exp ( - pow(( 2*local_y - 1-disp),power)  / (E) ) / (2*3.14159265359 * sqrt(E)); 
                // Co Variance
                C[j*N_U + i] *= sig_r1 * sig_r2 ;
            }

            else if(TDatabase::ParamDB->stddev_switch == 2){

              double sig_r1 = sin(-1.0*Pi*(2*actual_x-2))*sin(-1.0*Pi*(2*actual_y-2));
              double sig_r2 = sin(-1.0*Pi*(2*local_x-2))*sin(-1.0*Pi*(2*local_y-2));
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


//






    //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//======================================================================
// SystemMatrix construction and solution
//======================================================================  
    // Disc type: GALERKIN 
    // Solver: AMG_SOLVE (or) GMG  (or) DIRECT
    // SystemMatrix_Mean = new TSystemTBE2D(Velocity_FeSpace, Velocity_Mean, sol, rhs, Disctype, DIRECT);
    SystemMatrix = new TSystemTBE2D(Velocity_FeSpace, Velocity, sol, rhs, Disctype, DIRECT);


    // 2 parameters are needed for assembling (u1_old, u2_old)
    // FeFct[0] = u1_mean;
    // FeFct[1] = u2_mean; 
    // fesp[0] = Velocity_FeSpace;

    FeFct[0] = u1;
    FeFct[1] = u2; 
    fesp[0] = Velocity_FeSpace;
    BEaux = new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2, TimeNSN_ParamFct2,
                            TimeNSN_FEValues2, fesp, FeFct, TimeNSFct2, TimeNSFEFctIndex2, 
                            TimeNSFEMultiIndex2, TimeNSN_Params2, TimeNSBeginParam2);  
  
     // aux for calculating the error
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {
      BEaux_error =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
                             TimeNSN_ParamFct2,
                             TimeNSN_FEValues2,
                             fesp, FeFct,
                             TimeNSFct2,
                             TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
                             TimeNSN_Params2, TimeNSBeginParam2);     
      }
    // initilize the system matrix with the functions defined in Example file
    // last argument is aux that is used to pass additional fe functions (eg. mesh velocity)   
    // SystemMatrix_Mean->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue, BEaux, BEaux_error);
    // SystemMatrix_Mean->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue, BEaux, BEaux_error);
    SystemMatrix->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue, BEaux, BEaux_error);




    //Looping for Realization 
    for (int RealNo=0;RealNo<N_Realisations;RealNo++){

      cout << "Real No" << RealNo << endl;
    // assemble M, A matrices and rhs 
    // SystemMatrix_Mean->Assemble(sol, rhs);
    for(i=0;i<N_U;i++){
      sol[mappingArray[i]] = RealizationVector[RealNo+N_Realisations*i];
      sol[N_U+mappingArray[i]]=0;
      
    }
    SystemMatrix->Assemble(sol, rhs);
  
    // system for Mode
    // Disc type: GALERKIN 
    // Solver: AMG_SOLVE (or) GMG  (or) DIRECT
    // SystemMatrix_Mode = new TSystemTBE_Mode2D(Velocity_FeSpace, Velocity_Mode, sol, rhs, Disctype, DIRECT, Velocity_Mean);

    // // initilize the system matrix with the functions defined in Example file
    // // 2 parameters are needed for assembling (u1_old, u2_old)
    // FeFct_Mode[0] = u1_mean;
    // FeFct_Mode[1] = u2_mean; 
    // FeFct_Mode[2] = Velocity_Mode->GetComponent(0);
    // FeFct_Mode[3] = Velocity_Mode->GetComponent(1); 

    // fesp_mode[0] = VelocityMode_FeSpace;    
    // fesp_mode[1] = Velocity_FeSpace;    
    // BE_Modeaux = new TAuxParam2D(TimeNSN_FESpaces_Mode, TimeNSN_Fct_Mode, TimeNSN_ParamFct_Mode,
    //                       TimeNSN_FEValues_Mode, fesp_mode, FeFct_Mode, TimeNSFct_Mode, TimeNSFEFctIndex_Mode, 
    //                       TimeNSFEMultiIndex_Mode, TimeNSN_Params_Mode, TimeNSBeginParam_Mode);  

    // // last argument is aux that is used to pass additional fe functions (eg. mesh velocity)    
    // SystemMatrix_Mode->Init(LinCoeffs_Mode, BoundCondition, U1BoundValue, U2BoundValue, BEaux, BEaux_error);


//======================================================================
// produce outout
//======================================================================
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;   
  std::string str = std::to_string(RealNo);
  std::string filename = "Realization_Nr_" + std::to_string(RealNo);
  VtkBaseName = const_cast<char *>(filename.c_str());

  std::string name = "Realization _Number_" + std::to_string(int(RealNo));
  mkdir(name.c_str(), 0777); 
  Output = new TOutput2D(2, 2, 1, 1, Domain);

  //  Output->AddFEVectFunct(Velocity_Mean);
  Output->AddFEVectFunct(Velocity);
   
    if(TDatabase::ParamDB->WRITE_VTK)
     {
      os.seekp(std::ios::beg);
       if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
      img++;
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
   memset(AllErrors, 0, 7.*SizeOfDouble);
   
   // time loop starts
   while(TDatabase::TimeDB->CURRENTTIME< end_time)
    {                                               // time cycle
     m++;
     TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
     for(l=0;l<N_SubSteps;l++) // sub steps of fractional step theta
      {
       SetTimeDiscParameters(1);

      if(m==1)
       {
        OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
        OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
        OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
        OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
       }
       
      tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;
   
      // OutPut(endl << "CURRENT TIME: ");
      // OutPut(TDatabase::TimeDB->CURRENTTIME << endl);          
       
      //copy sol, rhs to olssol, oldrhs
      memcpy(oldrhs, rhs, N_Total_MeanDOF*SizeOfDouble);        
    
      // assemble only rhs, nonlinear matrix for NSE will be assemble in fixed point iteration
      // not needed if rhs is not time-dependent
      if(m!=1)
       { 
        //  SystemMatrix_Mean->AssembleA();
         SystemMatrix->AssembleA(); 
        }
      else
       { 
        //  SystemMatrix_Mean->Assemble(sol, rhs); 
         SystemMatrix->Assemble(sol, rhs);
          }

       
      //scale B matices and assemble NSE-rhs based on the \theta time stepping scheme 
      // SystemMatrix_Mean->AssembleSystMat(oldrhs, rhs, sol); 
      SystemMatrix->AssembleSystMat(oldrhs, rhs, sol);
      oldtau = tau;
         
      // calculate the residual
      // SystemMatrix_Mean->GetTBEResidual(sol, defect);
      SystemMatrix->GetTBEResidual(sol, defect);
      

    
      residual =  Ddot(N_Total_MeanDOF, defect, defect);
      // OutPut("Nonlinear iteration step   0");
      // OutPut(setw(14) << sqrt(residual) << endl);       

//====================================================================== 
//Solve the system
//Nonlinear iteration of fixed point type
//======================================================================
     for(j=1;j<=Max_It;j++)
      {   
       // Solve the NSE system
      //  SystemMatrix_Mean->Solve(sol);
      SystemMatrix->Solve(sol);
    
       // restore the mass matrix for the next nonlinear iteration      
      //  SystemMatrix_Mean->RestoreMassMat();     
      SystemMatrix->RestoreMassMat();     
    
       // assemble the system matrix with given aux, sol and rhs 
      //  SystemMatrix_Mean->AssembleANonLinear(sol, rhs);    
      SystemMatrix->AssembleANonLinear(sol, rhs);
          
       // assemble system mat, S = M + dt\theta_1*A
      //  SystemMatrix_Mean->AssembleSystMatNonLinear();        
      SystemMatrix->AssembleSystMatNonLinear();

       // get the residual
       memset(defect, 0, N_TotalDOF*SizeOfDouble);
      //  SystemMatrix_Mean->GetTBEResidual(sol, defect);  
      SystemMatrix->GetTBEResidual(sol, defect);        
 
    
       residual =  Ddot(N_Total_MeanDOF, defect, defect);
      //  OutPut("nonlinear iteration step " << setw(3) << j);
      //  OutPut(setw(14) << sqrt(residual) << endl); 
	
       if(sqrt(residual)<=limit)
       break;

      } // for(j=1;j<=Max_It;j++)

      
 
      // restore the mass matrix for the next time step          
      // SystemMatrix_Mean->RestoreMassMat();     
      SystemMatrix->RestoreMassMat();  
       
   } // for(l=0;l<N_SubSteps;
//====================================================================== 
// measure errors to known solution
//======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {          
      // SystemMatrix_Mean->MeasureErrors(ExactU1, ExactU2, AllErrors);
      SystemMatrix->MeasureErrors(ExactU1, ExactU2, AllErrors);

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
     if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)      
     if(TDatabase::ParamDB->WRITE_VTK)
      {
       os.seekp(std::ios::beg);
        if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());
       img++;
      }            
                
    } // while(TDatabase::TimeDB->CURRENTTIME< e

//======================================================================
// produce final outout
//======================================================================
    if(TDatabase::ParamDB->WRITE_VTK)
     {
      os.seekp(std::ios::beg);
       if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
      img++;
     }     
      TDatabase::TimeDB->CURRENTTIME = 0;
    } // Realization Loop
  CloseFiles();
  
  return 0;
} // end main
