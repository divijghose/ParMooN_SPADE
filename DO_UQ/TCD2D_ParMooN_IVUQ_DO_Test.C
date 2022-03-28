/**
 * @file TCD2D_ParMooN_IVUWQ_DO_Test.C
 * @brief Purpose:     Main program for scalar equations with new kernels of ParMooN. 
                       Features included in this main program - 
                       1. Monte Carlo Realization Generation for scalar quantity of interest
                       2. Initialization of Mean, Modes and Coefficients for solution of Dynamically Orthogonal System of Equations
                       3. 

 * @author Sashikumaar Ganesan
 * @author Divij Ghose
 * @author Thivin Anandh

 *  
 */

// =======================================================================
//
// Purpose:     Main program for scalar equations with new kernels of ParMooN. 
//              Features included in this main program - 
//              1. Monte Carlo Realization Generation for scalar quantity of interest
//              2. Initialization of Mean, Modes and Coefficients for solution of Dynamically Orthogonal System of Equations
//              3. <fill>
//
// Authors:      Sashikumaar Ganesan, Divij Ghose, Thivin Anandh
//
// History:     Implementation started on 10.03.2022

// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <SystemTCD2D.h>
#include <SystemTCD2D_ALE.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>
#include <TimeDiscRout.h>
#include<fstream>
#include<string> 
#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <mkl.h>
#include<cmath>
#include<random>
#include <FEVectFunct2D.h>



// =======================================================================
// include current example
// =======================================================================
// #include "../Examples/TCD_2D/exp.h"
#include "../Examples/DO_UQ/linear_advection_do_test.h"
// #include "../Examples/TCD_2D/SinCos1.h"
// #include "../Examples_All/TCD_2D/Time3.h"
// #include "../Examples/TCD_2D/exp_0.h"
//    #include "../Examples/TCD_2D/exp_2.h"
// #include "../Examples_All/TCD_2D/exp_1.h"
// #include "../Main_Users/Sashi/TCD_2D/Hemker.h"




int main(int argc, char *argv[])
{
    int i, j, l, m, N_SubSteps, ORDER, N_Cells, N_DOF, img = 1, N_G;
    int N_Active;

    double *sol, * rhs, *oldrhs, t1, t2, errors[5], Linfty;
    double tau, end_time, *defect, olderror, olderror1, hmin, hmax;

    double *solMean, *rhsMean, *old_rhsMean;

    bool UpdateStiffnessMat, UpdateRhs, ConvectionFirstTime;
    char *VtkBaseName;
    char *VtkBaseNameMean;
    char *VtkBaseNameMode;

    const char vtkdir[] = "VTK";

    TDomain *Domain;
    TDatabase *Database = new TDatabase();
    TFEDatabase2D *FEDatabase = new TFEDatabase2D();
    TCollection *coll;
    TFESpace2D *Scalar_FeSpace, *fesp[1];
    TFEFunction2D *Scalar_FeFunction_Mean;
    TFEFunction2D *Scalar_FeFunction;
    TOutput2D *Output;
    TOutput2D *OutputMean;
    TOutput2D *OutputMode;
    TSystemTCD2D *SystemMatrix;
    TAuxParam2D *aux;
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};

    std::ostringstream os;
    os << " ";//rename as os_mean?
    std::ostringstream os_mode;
	os_mode << " ";
    // ======================================================================
    // set the database values and generate mesh
    // ======================================================================
    // set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based
    Domain = new TDomain(argv[1]);

    if (TDatabase::ParamDB->PROBLEM_TYPE == 0)
        TDatabase::ParamDB->PROBLEM_TYPE = 2;
    OpenFiles();

    Database->WriteParamDB(argv[0]);
    Database->WriteTimeDB();
    ExampleFile();

    /* include the mesh from a mesh generator, for a standard mesh use the
   * build-in function. The GEOFILE describes the boundary of the domain. */
    if (TDatabase::ParamDB->MESH_TYPE == 0)
    {
        Domain->ReadGeo(TDatabase::ParamDB->GEOFILE);
        OutPut("PRM-GEO used for meshing !!!" << endl);
    } // ParMooN  build-in Geo mesh
    else if (TDatabase::ParamDB->MESH_TYPE == 1)
    {
        Domain->GmshGen(TDatabase::ParamDB->GEOFILE);
        OutPut("GMSH used for meshing !!!" << endl);
    }                                            //gmsh mesh
    else if (TDatabase::ParamDB->MESH_TYPE == 2) // ngle mesh
    {
        OutPut("Triangle.h used for meshing !!!" << endl);
        // TriaReMeshGen(Domain);
    }
    else
    {
        OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
        exit(0);
    }

#if defined(__HEMKER__) || defined(__BEAM__)
    TriaReMeshGen(Domain);
    TDatabase::ParamDB->UNIFORM_STEPS = 0;
#endif

    // refine grid up to the coarsest level
    for (i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
        Domain->RegRefineAll();

    // write grid into an Postscript file
    if (TDatabase::ParamDB->WRITE_PS)
        Domain->PS("Domain.ps", It_Finest, 0);

    // create output directory, if not already existing
    mkdir(vtkdir, 0777);
    const char modedir[] = "Modes";
    const char meandir[] = "Mean";
    const char coeffdir[] = "Coefficients";
    const char mcdir[] = "MonteCarlo";

    mkdir(meandir, 0777);
    mkdir(modedir, 0777);
    mkdir(coeffdir, 0777);
    mkdir(mcdir, 0777);

    //=========================================================================
    // construct all finite element spaces
    //=========================================================================
    ORDER = TDatabase::ParamDB->ANSATZ_ORDER;

    coll = Domain->GetCollection(It_Finest, 0);
    N_Cells = coll->GetN_Cells();
    OutPut("N_Cells (space) : " << N_Cells << endl);

    // fespaces for scalar equation
    Scalar_FeSpace = new TFESpace2D(coll, (char *)"fe space", (char *)"solution space",
                                    BoundCondition, 1, NULL);


    N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
    N_Active = Scalar_FeSpace->GetActiveBound();
    OutPut("dof all      : " << setw(10) << N_DOF << endl);
    OutPut("dof active   : " << setw(10) << N_Active << endl);


    ///////////////////////////////////////////////////////////////////////////////////////////////
    ////////// -------- REALISATION DATA GENERATION ----------------------------------------- //////
    ///////////////////////////////////////////////////////////////////////////////////////////////

    int N_Realisations  = TDatabase::ParamDB->REALIZATIONS;
    double LengthScale  = TDatabase::ParamDB->LENGTHSCALE;
    double EigenPercent = TDatabase::ParamDB->EIGENPERCENT;



    double *org_x_coord     = new double[N_DOF];
    double *org_y_coord     = new double[N_DOF];
    double *x_coord         = new double[N_DOF];
    double *y_coord         = new double[N_DOF];
    int *mappingArray       = new int[N_DOF];


    i=0;
    int N = pow(2,TDatabase::ParamDB->UNIFORM_STEPS  ) + 1;
    for ( int i = 0 ; i < N_DOF; i++)
    {
        int local_i = i/N;
        int local_j = i%N;
       
        x_coord[i] =  double(1.0/(N-1)) * local_i;
        y_coord[i] =  double(1.0/(N-1)) * local_j;
    }


    cout << " End File Read" <<endl;

    Scalar_FeSpace->GetDOFPosition(org_x_coord,org_y_coord);

    for ( int i=0 ; i < N_DOF; i++)   // Generated Values
    {  
        // get the generated Value
        double xx = x_coord[i];
        double yy = y_coord[i];
        bool foundFlag = false;

        for ( int j=0 ; j<N_DOF;j++)  // Actual parmooN Co-ordinates
        {  
            if(abs(xx - org_x_coord[j]) < 1e-10 &&  abs(yy - org_y_coord[j]) < 1e-10 )
            {
                mappingArray[i] = j;
                foundFlag = true;
            }
        }

        if(!foundFlag) cerr<< " DOF NOT FOUND FOR " << i << " position : " << setw(8) << org_x_coord[i]<<setw(8) <<org_y_coord[i] <<endl;
    }
     

    // int N_DOF =  N * N;
    double* x  =  new double[N_DOF];
    double* y  =  new double[N_DOF];

    for ( int i = 0 ; i < N_DOF; i++ )
    {
        int local_i = i/N;
        int local_j = i%N;

        x[i] =  double(1.0/(N-1)) * local_j;
        y[i] =  double(1.0/(N-1)) * local_i;
    }

   
    double *C = new double[N_DOF*N_DOF];  //MATRIX
    double *C1 = new double[N_DOF*N_DOF];  //MATRIX  - Corelation Matrix

    double norm = 0;
    for( int i =0  ; i < N_DOF ; i++ )
    {
        double actual_x = x[i];
        double actual_y = y[i];

        for ( int j=0 ; j < N_DOF ; j++)
        {
            double local_x = x[j];
            double local_y = y[j];

            double r = sqrt( pow((actual_x - local_x),2 ) + pow((actual_y - local_y),2 ));
            
            // CO -Relation
            C[j*N_DOF + i] = exp ( (- 1.0 * r )/ (LengthScale) );
            C1[j*N_DOF + i] = exp ( (- 1.0 * r )/ (LengthScale) );


            if(TDatabase::ParamDB->stddev_switch == 0)
            {
                double sig_r1 = exp (-1.0/(1.0 - pow(( 2*actual_x - 1),4) ) )  * exp ( -1.0/ ( 1 - pow(( 2*actual_y - 1),4) ) ) ;
                double sig_r2 = exp (-1.0/(1.0 - pow(( 2*local_x - 1),4) ) )  * exp ( -1.0/ ( 1 - pow(( 2*local_y - 1),4) ) ) ; 
            
                // Co Variance
                C[j*N_DOF + i] *= sig_r1 * sig_r2 * 5.0;
            }

            else if(TDatabase::ParamDB->stddev_switch == 1)
            {
                double E = TDatabase::ParamDB->stddev_denom;
                double disp = TDatabase::ParamDB->stddev_disp;
                double power = TDatabase::ParamDB->stddev_power;
                double sig_r1 = exp ( - pow( ( 2*actual_x - 1 - disp),power)  / (E) )  / (2*3.14159265359 * sqrt(E))  * exp ( - pow(( 2*actual_x - 1-disp),power)  / (E) )  / (2*3.14159265359 * sqrt(E)) ;
                double sig_r2 = exp ( - pow(( 2*local_x - 1 -disp),power)  / (E) )  / (2*3.14159265359 * sqrt(E))  * exp ( - pow(( 2*local_y - 1-disp),power)  / (E) ) / (2*3.14159265359 * sqrt(E)); 
                // Co Variance
                C[j*N_DOF + i] *= sig_r1 * sig_r2 ;
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

    for ( int i=0 ; i < N_DOF ; i++)
    {
        for ( int j=0 ; j < N_DOF ; j++)
        {
            fileo << C1[i*N_DOF + j] ;
            if(j!= N_DOF-1 ) fileo<<",";
        }
        fileo<<endl;
    }

    fileo.close();


    std::ofstream fileo_r;
    fileo.open("Covarriance.txt");

    for ( int i=0 ; i < N_DOF ; i++)
    {
        for ( int j=0 ; j < N_DOF ; j++)
        {
            fileo_r << C[i*N_DOF + j] ;
            if(j!= N_DOF-1 ) fileo_r<<",";
        }
        fileo_r<<endl;
    }

    fileo_r.close();

    // exit(0);
    ////////////////////////////////////////////////////// SVD ////////////////////////////////////////////
    // Declare SVD parameters
    MKL_INT m1 = N_DOF, n = N_DOF, lda = N_DOF, ldu = N_DOF, ldvt = N_DOF, info;
    double superb[std::min(N_DOF,N_DOF)-1];

    double* S = new double[N_DOF];
    double* U = new double[N_DOF*N_DOF];
    double* Vt = new double[N_DOF*N_DOF];
    cout << " REALISATIONS COMPUTED " <<endl;
    info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'A', 'A', m1, n, C, lda,
                        S, U, ldu, Vt, ldvt, superb );

    cout << endl <<endl;

    if( info > 0 ) {
                printf( "The algorithm computing SVD failed to converge.\n" );
                exit( 1 );
    }
   cout << " REALISATIONS COMPUTED " <<endl;
    int energyVal = 0;
    int temp = 0;


    double sumSingularVal = 0;
    for( int i=0;i<N_DOF;i++) sumSingularVal += S[i];

    double val = 0;
    for( energyVal =0 ; energyVal< N_DOF; energyVal++)
    {
        val += S[energyVal];
        temp++;
        if(val/sumSingularVal > 0.99) break;
    }

    cout << " MODES : "  << temp+1 <<endl;
    
    int modDim = temp+1;
       
    double* Ut = new double[N_DOF*modDim]();
    double* Z  = new double[N_Realisations*modDim]();

    double* RealizationVector = new double[N_DOF * N_Realisations]();
    double* RealizationVectorTemp = new double[N_DOF * N_Realisations]();


    // -------------- Generate Random Number Based on Normal Distribution -------------------------//
    int k=0;
    int skip = N_DOF - modDim;
    int count =0;
    for ( int i = 0 ; i < N_DOF*N_DOF ; i++ )
    {  
        // cout << "i val " << i <<endl;
        if(count < modDim )
        {
            Ut[k] =  U[i];

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


    for( int k = 0 ; k < modDim ; k++)
    {  
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<> d{0,1};

        double* norm1 = new double[N_Realisations];

        for(int n=0; n<N_Realisations; ++n) {
            Z[k*N_Realisations + n] =  S[k] * d(gen);
        }
    }

    cout << " N_Realisations : " << N_Realisations <<endl;
    cout << " MULT START "<<endl;
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,N_DOF,N_Realisations, modDim , 1.0, Ut,modDim,Z,N_Realisations,0.0,RealizationVector,N_Realisations);
    cout << " MULT DONE "<<endl;
    // printMatrix(SolutionVector, N_DOF,N_Realisations);

    // mkl_dimatcopy('R','T', N_DOF,N_Realisations,1.0,SolutionVector,N_DOF,N_Realisations);
    // cout << " COPY DONE "<<endl;

    cout << " REALISATIONS COMPUTED " <<endl;

    for(int i=0;i<N_DOF;i++){
        for(int j=0;j<N_Realisations;j++){
            RealizationVectorTemp[mappingArray[i]*N_Realisations+j]=RealizationVector[j+N_Realisations*i];
        }
    }

  memcpy(RealizationVector,RealizationVectorTemp,N_DOF*N_Realisations*SizeOfDouble);

    /////////////////////////////////////// -------- END OF REALISATION DATA SETS ------------ ////////////////////////////////////////////////////////////////

    ////////////////////////////////////// -------- START OF DO INITIALIZATION ------------ ////////////////////////////////////////////////////////////////

    double* MeanVector = new double[N_DOF * 1](); //overline{C}_{dof} = \sum_{i=1}^{N_Realisations}(C^{i}_{dof})/N_Realisations
    for(int i=0; i<N_DOF; ++i) {
    for(int j = 0; j < N_Realisations; ++j){
                MeanVector[i] +=  (RealizationVector[i*N_Realisations+j]/N_Realisations);
            }
    }

    double* PerturbationVector = new double[N_DOF * N_Realisations](); // \hat{C}^{i}_{dof} = C^{i}_{dof} - \overline{C}_{dof}
    for(int i = 0; i < N_DOF; ++i){
    for(int j = 0; j < N_Realisations; ++j){
        PerturbationVector[i*N_Realisations+j] = RealizationVector[i*N_Realisations+j] - MeanVector[i];
    }
    }

//================================================================================================
/////////////////////////////DO - Initialization SVD//////////////////////////////////////////////
//================================================================================================
// Declare SVD parameters
int minDim = std::min(N_DOF,N_Realisations);
MKL_INT mDO = N_DOF, nDO = N_Realisations, ldaDO = N_Realisations, lduDO = minDim, ldvtDO = N_Realisations, infoDO;
double superbDO[minDim-1];

double* PerturbationVectorCopy = new double[N_DOF * N_Realisations]();
memcpy(PerturbationVectorCopy,PerturbationVector,N_DOF*N_Realisations*SizeOfDouble);

double* Sg = new double[minDim];
double* L = new double[N_DOF*minDim];
double* Rt = new double[minDim*N_Realisations];

infoDO = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'S', 'N', mDO, nDO, PerturbationVectorCopy, ldaDO,
					Sg, L, lduDO, Rt, ldvtDO, superbDO );


if( infoDO > 0 ) {
	printf( "The algorithm computing SVD for DO failed to converge.\n" );
	exit( 1 );
}
cout << " DO SVD COMPUTED " <<endl;

//////////////////////////////////////////// DO - SVD End/////////////////////////////// 

///////DO - Subspace dimension calculation /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
double SVPercent = TDatabase::ParamDB->SVPERCENT;
int s = 0;
double valDO = 0.0;
double sumSingularValDO = 0.0;
for( int i=0;i<minDim;i++) {
    sumSingularValDO += Sg[i];
}
while( valDO/sumSingularValDO < SVPercent)
{
    valDO += Sg[s];
    s++;
}

cout << " SUBSPACE DIMENSION : "  << s+1 <<endl;

int subDim = s+1;
////////Subspace dimension calculated//////////////////

/////Projection Matrix///////////
////////////////////////////////
cout << " Min DIMENSION : "  << minDim <<endl;
double* ProjectionVector = new double[N_Realisations * minDim]();
cout << "PROJ VECTOR MULT START "<<endl;
cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,N_Realisations,minDim,N_DOF,1.0,PerturbationVector,N_Realisations,L,minDim,0.0,ProjectionVector,minDim);
cout << "PROJ VECTOR MULT DONE "<< endl;

/// Initialize Coefficient Matrix - First subDim columns of Projection Matrix ////////////////////
double* CoeffVector = new double[N_Realisations * subDim]();
// memcpy(CoeffVector, ProjectionVector, N_Realisations*subDim*SizeOfDouble); //For ColMajor storage -wrong!!
// for (int i=0;i<N_Realisations;i++){
// 	for (int j=0;j<subDim;j++){
// 		CoeffVector[i*subDim+j] = ProjectionVector[i*minDim+j];
// 	}
// }
for(int i=0;i<N_Realisations;i++){
    for(int j=0;j<subDim;j++){
        CoeffVector[j*N_Realisations+i] = ProjectionVector[i*minDim+j]; // CoeffVector in Col Major form 
    }
}

////////////Initialize Mode Vector - First subDim columns of Left Singular Vector//////////////////
double* ModeVector = new double[N_DOF* subDim]();
// memcpy(ModeVector, L, N_DOF*subDim*SizeOfDouble);//For ColMajor storage
// for (int i=0;i<N_DOF;i++){
// 	for (int j=0;j<subDim;j++){
// 		ModeVector[i*subDim+j] = ProjectionVector[i*minDim+j];
// 	}
// }
for(int i=0;i<N_DOF;i++){
    for(int j=0;j<subDim;j++){
        ModeVector[j*N_DOF+i] = L[i*minDim+j]; // ModeVector in Col Major form 
    }
}

////////////////////////////////////////////DO - Initialization Ends//////////////////////////////////////
///////================================================================================//////////////////






//======================================================================
// construct all finite element functions
//======================================================================
sol = new double[N_DOF];
rhs = new double[N_DOF];
oldrhs = new double[N_DOF];

memset(sol, 0, N_DOF * SizeOfDouble);
memset(rhs, 0, N_DOF * SizeOfDouble);
memset(oldrhs, 0, N_DOF * SizeOfDouble);

solMean = new double[N_DOF];
rhsMean = new double[N_DOF];
old_rhsMean = new double[N_DOF];

memset(solMean, 0, N_DOF * SizeOfDouble);
memset(rhsMean, 0, N_DOF * SizeOfDouble);
memset(old_rhsMean, 0, N_DOF * SizeOfDouble);

Scalar_FeFunction_Mean = new TFEFunction2D(Scalar_FeSpace, (char *)"C_Mean", (char *)"sol", solMean, N_DOF);

Scalar_FeFunction_Mean->Interpolate(InitialCondition);
for(int i=0;i<N_DOF;i++){
    // solMean[mappingArray[i]]=MeanVector[i];
    solMean[i]=MeanVector[i];
}

// TDatabase::ParamDB->COVARIANCE_MATRIX_DO = CalcCovarianceMatx(CoeffVector,N_Realisations,subDim);//Not needed for linear advection

	//======================================================================
	// /DO - SystemMatrix construction and solution
	//======================================================================
	// Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
	// Solver: AMG_SOLVE (or) GMG  (or) DIRECT
	TSystemTCD2D* SystemMatrix_Mean = new TSystemTCD2D(Scalar_FeSpace, GALERKIN, DIRECT);
    TSystemTCD2D* SystemMatrix_Mode =  new TSystemTCD2D(Scalar_FeSpace, GALERKIN, DIRECT);

	// initilize the system matrix with the functions defined in Example file
	SystemMatrix_Mean->Init(DO_Mean_Equation_Coefficients, BoundCondition, BoundValue);
	SystemMatrix_Mode->Init(DO_Mode_Equation_Coefficients, BoundCondition, BoundValue);

    // 

	   // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //------------------------------------------ MEAN EQUATION SETUP -----------------------------------------------------//
    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//

    int N_terms_Mean = 3;  // Number of Shape function derivatives required ( in this case 3, N, NX, NY )
    MultiIndex2D Derivatives_MatrixARhs_Mean[3] = { D00, D10, D01};
    int SpacesNumbers_MatrixARhs_Mean[3] = { 0, 0, 0  };
    int N_Matrices_MatrixARhs_Mean = 1;
    int RowSpace_MatrixARhs_Mean[1] = { 0 };
    int ColumnSpace_MatrixARhs_Mean[1] = { 0 };
    int N_Rhs_MatrixARhs_Mean = 1;
    int RhsSpace_MatrixARhs_Mean[1] = { 0 };


    SystemMatrix_Mean->Init_WithDiscreteform(DO_Mean_Equation_Coefficients,BoundCondition, BoundValue,"DO_LINEAR_MEAN", "DO_LINEAR_MEAN",
                                                N_terms_Mean, Derivatives_MatrixARhs_Mean, SpacesNumbers_MatrixARhs_Mean,
                                                N_Matrices_MatrixARhs_Mean, N_Rhs_MatrixARhs_Mean,
                                                RowSpace_MatrixARhs_Mean,ColumnSpace_MatrixARhs_Mean, RhsSpace_MatrixARhs_Mean,
                                                DO_Mean_Equation_Assembly, DO_Mean_Equation_Coefficients,
											NULL);

    // Aux Setup for the RHS -- There is no Aux for the Mean equation, So set the values as NULL
    fesp[0] = Scalar_FeSpace;
    aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

    /* -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-
    --------------------------------------[[[ END  ]]] MEAN EQUATION SETUP -----------------------------------------------------
     -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0- */
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	 // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //------------------------------------------ CO EFFICIENT EQUATION SETUP -----------------------------------------------------//
    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//

	TFEVectFunct2D* FeVector_Coefficient = new TFEVectFunct2D(Scalar_FeSpace, (char *)"sol", (char *)"sol",CoeffVector, N_Realisations, subDim);

		 // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //------------------------------------------ CO EFFICIENT EQUATION SETUP END -------------------------------------------------------//
    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//


	// Aux Setup for the RHS -- There is no Aux for the Mean equation, So set the values as NULL
    fesp[0] = Scalar_FeSpace;
    aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	 // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //-------------------------------------- MODE EQUATION SETUP -----------------------------------------------------//
    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
	 int N_terms_Mode = 3;  // Number of Shape function derivatives required ( in this case 3, N, NX, NY )
    MultiIndex2D Derivatives_MatrixARhs_Mode[3] = { D00, D10, D01};
    int SpacesNumbers_MatrixARhs_Mode[3] = { 0, 0, 0  };
    int N_Matrices_MatrixARhs_Mode = 1;
    int RowSpace_MatrixARhs_Mode[1] = { 0 };
    int ColumnSpace_MatrixARhs_Mode[1] = { 0 };
    int N_Rhs_MatrixARhs_Mode = 1;
    int RhsSpace_MatrixARhs_Mode[1] = { 0 };
	SystemMatrix_Mode->Init_WithDiscreteform(DO_Mode_Equation_Coefficients,BoundCondition, BoundValue,"DO_LINEAR_Mode", "DO_LINEAR_Mode",
                                                N_terms_Mode, Derivatives_MatrixARhs_Mode, SpacesNumbers_MatrixARhs_Mode,
                                                N_Matrices_MatrixARhs_Mode, N_Rhs_MatrixARhs_Mode,
                                                RowSpace_MatrixARhs_Mode,ColumnSpace_MatrixARhs_Mode, RhsSpace_MatrixARhs_Mode,
                                                DO_Mode_Equation_Assembly, DO_Mode_Equation_Coefficients,
                                                NULL);


	double* solMode = new double[N_DOF* subDim]();
	double* rhsMode = new double[N_DOF* subDim]();
	double* old_rhsMode = new double[N_DOF]();

    for(int j=0;j<subDim;j++){
        for(int i =0; i<N_DOF;i++){
            // solMode[j*N_DOF+mappingArray[i]] = ModeVector[j*N_DOF+i];
            solMode[j*N_DOF+i] = ModeVector[j*N_DOF+i];

        }
    }




// 	// double* ModeVector_OldRHS = new double[N_DOF]();
	TFEVectFunct2D* FEFVector_Mode = new TFEVectFunct2D(Scalar_FeSpace, (char *)"C_Mode", (char *)"sol",solMode, N_DOF, subDim);

// // Set up a FE VECT FUNCTION TO STORE ALL THE Components of CTilde  

//     // TFEVectFunct2D* linModesFeVectFunct = 

    int TimeLinear_FESpaces_DO = 1;
	int TimeLinear_Fct_DO = 1; // \tilde(C)
	int TimeLinear_ParamFct_DO = 1;
	int TimeLinear_FEValues_DO = 3;
	int TimeLinear_Params_DO = 3;
	int TimeNSFEFctIndex_DO[3] = {0, 0, 0};  
	MultiIndex2D TimeNSFEMultiIndex_DO[3] = {D00, D01, D10};
	ParamFct *TimeNSFct_DO[1] 		= { DO_Mode_RHS_Aux_Param };
	int TimeNSBeginParam_DO[1] 		= { 0 };

    TFEFunction2D* fefct_RHS[4];
	TFESpace2D* fesp_RHS[2];

// ;
//     // Set Up Aux  Param for the given MOde 
//     // The mode equation needs the entire values of the C_tilde matrix array 
//     // TAuxParam2D* aux_RHS_DO = new TAuxParam2D (	TimeLinear_FESpaces_DO, <FE VECT FUNCTION>, TimeLinear_ParamFct_DO,
// 	// 											TimeLinear_FEValues_DO,
// 	// 											fesp_RHS,
// 	// 											TimeNSFct_DO,
// 	// 											TimeNSFEMultiIndex_DO,
// 	// 											TimeLinear_Params_DO, TimeNSBeginParam_DO);

	TAuxParam2D* aux_RHS_DO = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
	  // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
    //--------------------------------------[[[ END  ]]] MODE EQUATION SETUP -----------------------------------------------------//
    // -0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0---0-0--0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0--00-0-0-0-0-0-0-0-0-0-0-0--0-0-0-0-//
// assemble the system matrix with given aux, sol and rhs
	// aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
	// otherwise, just pass with NULL
	SystemMatrix_Mean->AssembleMRhs(NULL, solMean, rhsMean);
	
	SystemMatrix_Mode->AssembleMRhs(NULL, solMode, rhsMode);




	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//======================================================================
	// produce outout at t=0
	//======================================================================
	// VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
    std::string strMean= std::to_string(N_Realisations);
    std::string filenameMean = "Mean_NRealisations_" + std::to_string(N_Realisations);
    VtkBaseNameMean = const_cast<char *>(filenameMean.c_str());

    std::string strMode= std::to_string(N_Realisations);
    std::string filenameMode = "Mode_NRealisations_" + std::to_string(N_Realisations);
    VtkBaseNameMode = const_cast<char *>(filenameMode.c_str());


    std::string fileoutCoeff; 
    std::string fileoutMean; 
    std::string fileoutMode; 
    std::string fileoutMC;    


	OutputMean = new TOutput2D(2, 2, 1, 1, Domain);
	OutputMean->AddFEFunction(Scalar_FeFunction_Mean);

    OutputMode = new TOutput2D(2, 2, 1, 1, Domain);
	OutputMode->AddFEVectFunct(FEFVector_Mode);


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
	
	
	 if(TDatabase::ParamDB->WRITE_VTK)
     {
      os.seekp(std::ios::beg);
       if(modeimg<10) os <<  "VTK/"<<VtkBaseNameMode<<".0000"<<modeimg<<".vtk" << ends;
         else if(modeimg<100) os <<  "VTK/"<<VtkBaseNameMode<<".000"<<modeimg<<".vtk" << ends;
          else if(modeimg<1000) os <<  "VTK/"<<VtkBaseNameMode<<".00"<<modeimg<<".vtk" << ends;
           else if(modeimg<10000) os <<  "VTK/"<<VtkBaseNameMode<<".0"<<modeimg<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseNameMode<<"."<<modeimg<<".vtk" << ends;
      OutputMode->WriteVtk(os.str().c_str());
      modeimg++;
     }     
	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	coll->GetHminHmax(&hmin, &hmax);
	OutPut("h_min : " << hmin << " h_max : " << hmax << endl);

	//======================================================================
	// time disc loop
	//======================================================================
	// parameters for time stepping scheme
	m = 0;
	N_SubSteps = GetN_SubSteps();
	end_time = TDatabase::TimeDB->ENDTIME;

	UpdateStiffnessMat = FALSE; //check BilinearCoeffs in example file
	UpdateRhs = FALSE;			//check BilinearCoeffs in example file
	ConvectionFirstTime = TRUE;
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


        //======================================================================
        // time disc loop
        //======================================================================
        // parameters for time stepping scheme
        m = 0;
        N_SubSteps = GetN_SubSteps();
        end_time = TDatabase::TimeDB->ENDTIME;

        UpdateStiffnessMat = TRUE; //check BilinearCoeffs in example file
        UpdateRhs = TRUE;           //check BilinearCoeffs in example file
        ConvectionFirstTime = TRUE;

        // time loop starts
        while (TDatabase::TimeDB->CURRENTTIME < end_time)
        {
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

                tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
                TDatabase::TimeDB->CURRENTTIME += tau;

                OutPut(endl<< "CURRENT TIME: ");
                OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

                // copy rhs to oldrhs
                memcpy(old_rhsMean, rhsMean, N_DOF * SizeOfDouble);

                // unless the stiffness matrix or rhs change in time, it is enough to
                // assemble only once at the begning
				SystemMatrix_Mean->AssembleARhs(NULL, solMean, rhsMean);
          
				for ( int subSpaceNum = 0 ; subSpaceNum < subDim; subSpaceNum++)
			{
                // solve the system matrix
                DO_CoEfficient(Scalar_FeSpace, FEFVector_Mode ,FeVector_Coefficient, subDim, subSpaceNum,N_Realisations);
				double* modeSolution_i = solMode + subSpaceNum*N_DOF; // This works for column major 
				double* modeSolution_rhs = rhsMode + subSpaceNum*N_DOF; // This works for column major 
				
				memcpy(old_rhsMode, modeSolution_rhs, N_DOF * SizeOfDouble);
				
				SystemMatrix_Mode->AssembleARhs(NULL, modeSolution_i, modeSolution_rhs);
				
				// get the UPDATED RHS VALUE FROM FUNCTION
				DO_Mode_RHS(Scalar_FeSpace,FEFVector_Mode, subDim, modeSolution_rhs,subSpaceNum);
				
				SystemMatrix_Mode->AssembleSystMat(old_rhsMode, modeSolution_i, modeSolution_rhs, modeSolution_i);
				SystemMatrix_Mode->Solve(modeSolution_i, modeSolution_rhs);
				SystemMatrix_Mode->RestoreMassMat();
                
                }//subSpaceNumLoop



        fileoutCoeff= "Coefficients/Coeff_NRealisations_" + std::to_string(N_Realisations)+"_t"+std::to_string(m);
        std::ofstream fileCoeff;
        fileCoeff.open(fileoutCoeff);

        for (int i = 0; i < N_Realisations; i++)
        {
            for (int j = 0; j < subDim; j++)
            {
                fileCoeff << CoeffVector[j*subDim + i];
                if (j != subDim - 1)
                    fileCoeff << ",";
            }
            fileCoeff << endl;
        }

        fileCoeff.close();

        fileoutMode= "Modes/Mode_NRealisations_" + std::to_string(N_Realisations)+"_t"+std::to_string(m);
        std::ofstream fileMode;
        fileMode.open(fileoutMode);

        for (int i = 0; i < N_DOF; i++)
        {
            for (int j = 0; j < subDim; j++)
            {
                fileMode << solMode[j*subDim + i];
                if (j != subDim - 1)
                    fileMode << ",";
            }
            fileMode << endl;
        }

        fileMode.close();


        fileoutMean= "Mean/Mean_NRealisations_" + std::to_string(N_Realisations)+"_t"+std::to_string(m);
        std::ofstream fileMean;
        fileMean.open(fileoutMean);

        for (int i = 0; i < N_DOF; i++)
        {
            
                fileMean << MeanVector[i]<<",";
                fileMean << endl;
        }

        fileMean.close();

				
                
            SystemMatrix_Mean->AssembleSystMat(old_rhsMean, solMean, rhsMean, solMean);
            ////  -- 
			SystemMatrix_Mean->Solve(solMean, rhsMean);
			SystemMatrix_Mean->RestoreMassMat(); 




                // restore the mass matrix for the next time step
                // unless the stiffness matrix or rhs change in time, it is not necessary to assemble the system matrix in every time step
             

            } // for(l=0;l< N_SubSteps;l++)
			//======================================================================
		// produce outout
		//======================================================================
		if (m == 1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
		{
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

			

			if(TDatabase::ParamDB->WRITE_VTK)
			{
				os.seekp(std::ios::beg);
				if(modeimg<10) os <<  "VTK/"<<VtkBaseNameMode<<".0000"<<modeimg<<".vtk" << ends;
					else if(modeimg<100) os <<  "VTK/"<<VtkBaseNameMode<<".000"<<modeimg<<".vtk" << ends;
					else if(modeimg<1000) os <<  "VTK/"<<VtkBaseNameMode<<".00"<<modeimg<<".vtk" << ends;
					else if(modeimg<10000) os <<  "VTK/"<<VtkBaseNameMode<<".0"<<modeimg<<".vtk" << ends;
						else  os <<  "VTK/"<<VtkBaseNameMode<<"."<<modeimg<<".vtk" << ends;
				OutputMode->WriteVtk(os.str().c_str());
				modeimg++;
			} 


		} //produce output loop end
            //======================================================================
		// measure errors to known solution
		//======================================================================
		if (TDatabase::ParamDB->MEASURE_ERRORS)
		{
            fesp[0] = Scalar_FeSpace;
            aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
			Scalar_FeFunction_Mean->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, DO_Mean_Equation_Coefficients, aux, 1, fesp, errors);

			OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
			OutPut(" L2: " << errors[0]);
			OutPut(" H1-semi: " << errors[1] << endl);

			errors[3] += (errors[0] * errors[0] + olderror * olderror) * TDatabase::TimeDB->TIMESTEPLENGTH / 2.0;
			olderror = errors[0];
			OutPut(TDatabase::TimeDB->CURRENTTIME << " L2(0,T;L2) " << sqrt(errors[3]) << " ");

			errors[4] += (errors[1] * errors[1] + olderror1 * olderror1) * TDatabase::TimeDB->TIMESTEPLENGTH / 2.0;
			OutPut("L2(0,T;H1) " << sqrt(errors[4]) << endl);
			olderror1 = errors[1];

			if (Linfty < errors[0])
				Linfty = errors[0];

			OutPut("Linfty " << Linfty << endl);
		} //  if(TDatabase::ParamDB->MEASURE_ERRORS)

        } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

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

			

			if(TDatabase::ParamDB->WRITE_VTK)
			{
				os.seekp(std::ios::beg);
				if(modeimg<10) os <<  "VTK/"<<VtkBaseNameMode<<".0000"<<modeimg<<".vtk" << ends;
					else if(modeimg<100) os <<  "VTK/"<<VtkBaseNameMode<<".000"<<modeimg<<".vtk" << ends;
					else if(modeimg<1000) os <<  "VTK/"<<VtkBaseNameMode<<".00"<<modeimg<<".vtk" << ends;
					else if(modeimg<10000) os <<  "VTK/"<<VtkBaseNameMode<<".0"<<modeimg<<".vtk" << ends;
						else  os <<  "VTK/"<<VtkBaseNameMode<<"."<<modeimg<<".vtk" << ends;
				OutputMode->WriteVtk(os.str().c_str());
				modeimg++;
			} 

    cout <<"Subspace Dimension = " << subDim << endl;
    CloseFiles();


//////////////////////// Monte Carlo Runs /////////////////////////////////////////////
    TDatabase::TimeDB->CURRENTTIME = 0;
    sol = new double[N_DOF];
    rhs = new double[N_DOF];
    oldrhs = new double[N_DOF];

    memset(sol, 0, N_DOF * SizeOfDouble);
    memset(rhs, 0, N_DOF * SizeOfDouble);

    Scalar_FeFunction = new TFEFunction2D(Scalar_FeSpace, (char *)"sol", (char *)"sol", sol, N_DOF);    
    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
    
    Output = new TOutput2D(2, 2, 1, 1, Domain);
    Output->AddFEFunction(Scalar_FeFunction);



    std::ofstream fileout;
    std::string name = "init" ;
    fileout.open(name);

    std::ofstream fileout_final;
    std::string name1 = "Final" ;
    fileout_final.open(name1);

    //======================================================================
    // SystemMatrix construction and solution
    //======================================================================
    // Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
    // Solver: AMG_SOLVE (or) GMG  (or) DIRECT
    SystemMatrix = new TSystemTCD2D(Scalar_FeSpace, GALERKIN, DIRECT);

    // initilize the system matrix with the functions defined in Example file
    SystemMatrix->Init(BilinearCoeffs, BoundCondition, BoundValue);



    // Setup array for random number
    srand(time(NULL));
    int N_samples  = 100;
    int* indexArray =  new int[N_samples];
    for ( int i =0 ; i < N_samples;i++)
        indexArray[i] = rand() % N_DOF;
    


    for ( int RealNo=0 ; RealNo < N_Realisations; RealNo++)
    {
        // cout << " ============================================================================================================= " <<endl;
        cout <<  " Real no " <<  RealNo  <<endl;
        // cout << " ============================================================================================================= " <<endl;

		
		
		
		
        // Scalar_FeFunction->Interpolate(InitialCondition);
        std::string str= std::to_string(RealNo);
        std::string filename = "Realization_Nr_" + std::to_string(RealNo);
		VtkBaseName = const_cast<char *>(filename.c_str());
        img=0;

        mkdir(filename.c_str(), 0777);

        for ( int i=0 ; i < N_DOF; i++)
            sol[i] = RealizationVector[RealNo  +  N_Realisations * i];
        
        os.seekp(std::ios::beg);
        if (img < 10)
            os << filename.c_str() << "/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
        else if (img < 100)
            os << filename.c_str() << "/" << VtkBaseName<< ".000" << img << ".vtk" << ends;
        else if (img < 1000)
            os << filename.c_str() << "/" << VtkBaseName<< ".00" << img << ".vtk" << ends;
        else if (img < 10000)
           os << filename.c_str() << "/" << VtkBaseName<< ".0" << img << ".vtk" << ends;
        else
            os << filename.c_str() << "/" << VtkBaseName << "." << img << ".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
        img++;

        for ( int k = 0 ; k < N_samples;k++)
        {
            fileout << sol[indexArray[k]];
            if(k != N_DOF - 1) fileout <<",";
        }
        fileout<<endl;

        // assemble the system matrix with given aux, sol and rhs
        // aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
        // otherwise, just pass with NULL
        SystemMatrix->AssembleMRhs(NULL, sol, rhs);

        //======================================================================
        // time disc loop
        //======================================================================
        // parameters for time stepping scheme
        m = 0;
        N_SubSteps = GetN_SubSteps();
        end_time = TDatabase::TimeDB->ENDTIME;

        UpdateStiffnessMat = TRUE; //check BilinearCoeffs in example file
        UpdateRhs = TRUE;           //check BilinearCoeffs in example file
        ConvectionFirstTime = TRUE;

        // time loop starts
        while (TDatabase::TimeDB->CURRENTTIME < end_time)
        {
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

                tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
                TDatabase::TimeDB->CURRENTTIME += tau;

                // OutPut(endl<< "CURRENT TIME: ");
                // OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

                //copy rhs to oldrhs
                memcpy(oldrhs, rhs, N_DOF * SizeOfDouble);

                // unless the stiffness matrix or rhs change in time, it is enough to
                // assemble only once at the begning
                if (UpdateStiffnessMat || UpdateRhs || ConvectionFirstTime)
                {
                    SystemMatrix->AssembleARhs(NULL, sol, rhs);

                    // M:= M + (tau*THETA1)*A
                    // rhs: =(tau*THETA4)*rhs +(tau*THETA3)*oldrhs +[M-(tau*THETA2)A]*oldsol
                    // note! sol contains only the previous time step value, so just pass
                    // sol for oldsol
                    SystemMatrix->AssembleSystMat(oldrhs, sol, rhs, sol);
                    ConvectionFirstTime = FALSE;
                }

                // solve the system matrix
                SystemMatrix->Solve(sol, rhs);

            fileoutMC= "MonteCarlo/MC_NRealisations_" + std::to_string(N_Realisations)+"_t"+std::to_string(m);
            std::ofstream fileMC;
            fileMC.open(fileoutMC,std::fstream::app);

        for (int i = 0; i < N_DOF; i++)
        {
            fileMC << sol[i] << ",";
        }
        fileMC << endl;
        cout << "Realization "<< RealNo <<" added" <<endl;

        fileMC.close();

                // restore the mass matrix for the next time step
                // unless the stiffness matrix or rhs change in time, it is not necessary to assemble the system matrix in every time step
                if (UpdateStiffnessMat || UpdateRhs)
                {
                    SystemMatrix->RestoreMassMat();
                }

                if (TDatabase::ParamDB->WRITE_VTK)
                {
                    os.seekp(std::ios::beg);
                    if (img < 10)
                        os << filename.c_str() << "/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
                    else if (img < 100)
                        os << filename.c_str() << "/" << VtkBaseName << ".000" << img << ".vtk" << ends;
                    else if (img < 1000)
                        os << filename.c_str() << "/" << VtkBaseName << ".00" << img << ".vtk" << ends;
                    else if (img < 10000)
                        os << filename.c_str() << "/" << VtkBaseName << ".0" << img << ".vtk" << ends;
                    else
                        os << filename.c_str() << "/" << VtkBaseName << "." << img << ".vtk" << ends;
                    Output->WriteVtk(os.str().c_str());
                    img++;
                }

            } // for(l=0;l< N_SubSteps;l++)

        } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

        //======================================================================
        // produce final outout
        //======================================================================

        os.seekp(std::ios::beg);
        if (img < 10)
            os << filename.c_str() << "/" << VtkBaseName << ".0000" << img << ".vtk" << ends;
        else if (img < 100)
            os << filename.c_str() << "/" << VtkBaseName << ".000" << img << ".vtk" << ends;
        else if (img < 1000)
            os << filename.c_str() << "/" << VtkBaseName << ".00" << img << ".vtk" << ends;
        else if (img < 10000)
            os << filename.c_str()<< "/" << VtkBaseName << ".0" << img << ".vtk" << ends;
        else
            os << filename.c_str() << "/" << VtkBaseName << "." << img << ".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
        img++;

        // cout << " Solution Norm After: " << Ddot(N_DOF,sol,sol) <<endl;

        for ( int k = 0 ; k < N_samples;k++)
        {
            fileout_final << sol[indexArray[k]];
            if(k != N_DOF - 1) fileout_final <<",";
        }
        fileout_final<<endl;
        

        // set Current Time as Zero
        TDatabase::TimeDB->CURRENTTIME = 0;
        // delete SystemMatrix;

    }
    fileout.close();

    CloseFiles();



/////////////////////// Monte Carlo Run Ends /////////////////////////////////////////

    return 0;
} // end main
