

/**
 * @file navierstokes_do_test.h
 * @brief Purpose:     Example file for solving the set of dynamically orthogonal field
					   equations for Navier-Stokes' equation.
					   Features included in this example file :
					   1. Definition of boundary conditions, boundary values,
					   initial condition and bilinear coefficients
					   2. Assembly functions for Mean Equation, Mode Equation
					   and Coefficient Equation
					   3. Assembly function for RHS of Mode Equation

 * @authors Sashikumaar Ganesan
 * @authors Divij Ghose
 * @authors Thivin Anandh
 * @bug No known bugs
 */

// ===========================================================================//
//  Dynamically Orthogonal Field Equation Solution of Navier-Stokes' Equation //
// ===========================================================================//

// =======================================================================
//
// Purpose:     Example file for solving the set of dynamically orthogonal field
//              equations for Navier-Stokes' equation.
//              Features included in this example file -
//              1. Definition of boundary conditions, boundary values,
//                 initial condition and bilinear coefficients
//              2. Assembly functions for Mean Equation, Mode Equation
//                 and Coefficient Equation
//              3. Assembly function for RHS of Mode Equation
//
// Authors:      Sashikumaar Ganesan, Thivin Anandh, Divij Ghose
//
// History:     1> First iteration implemented on 30.03.2022
//		          2> Bug fixes on 24.03.2022

// =======================================================================

void ExampleFile()
{

	OutPut("Example: Navier Stokes' Equation - Dynamically Orthogonal Field Equation Solution" << endl);
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1Mean(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialU2Mean(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialPMean(double x, double y, double *values)
{
  values[0] = 0;
}


void InitialU1(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = 0.0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
	cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0: 
            value=0;
            break;
    case 1: 
            value=0;
            break;
    case 2: 
            if(abs(Param - 0) < 1e-6 || abs(Param - 1.0) < 1e-6 )
            value=0; // top moving side velocity
            else 
              value = 1.0;
            break;
    case 3: 
            value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

void BoundConditionMean(int i, double t, BoundCond &cond)
{
	cond = DIRICHLET;
}

void U1BoundValueMean(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0: 
            value=0;
            break;
    case 1: 
            value=0;
            break;
    case 2: 
            if(abs(Param - 0) < 1e-6 || abs(Param - 1.0) < 1e-6 )
            value=0; // top moving side velocity
            else 
              value = 1.0;
            break;
    case 3: 
            value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U2BoundValueMean(int BdComp, double Param, double &value)
{
  value = 0;
}

void BoundConditionMode(int i, double t, BoundCond &cond)
{
	cond = DIRICHLET;
}

void U1BoundValueMode(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0: 
            value=0;
            break;
    case 1: 
            value=0;
            break;
    case 2: 
            if(abs(Param - 0) < 1e-6 || abs(Param - 1.0) < 1e-6 )
            value=0; // top moving side velocity
            else 
              value = 1.0;
            break;
    case 3: 
            value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U2BoundValueMode(int BdComp, double Param, double &value)
{
  value = 0;
}



// data on each quadrature point for mean
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  int i;
  double *coeff, x, y; 
  static double eps=double(1.0/TDatabase::ParamDB->RE_NR);

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = eps;

    coeff[1] = 0;  // f1
    coeff[2] = 0;  // f2
  }
}


void DO_Mean_Equation_Coefficients(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  int i;
  double *coeff, x, y; 
  static double eps=double(1.0/TDatabase::ParamDB->RE_NR);

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = eps;

    coeff[1] = 0;  // f1
    coeff[2] = 0;  // f2
  }
}

void DO_Mode_Equation_Coefficients(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  int i;
  double *coeff, x, y; 
  static double eps=double(1.0/TDatabase::ParamDB->RE_NR);

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = eps;

    coeff[1] = 0;  // f1
    coeff[2] = 0;  // f2
  }
}



/**
 * @brief Routine to calculate the covariance matrix of coefficients and store it in TDatabas::ParamDB->COVARIANCE_MATRIX_DO
 * \f{equation}{
  C = \frac{1}{N_{R}-1}\cdot\phi^{T}\phi
\f}
 *
 * @param Vector Coefficient Matrix
 * @param height Number of Realizations, N_R (number of rows of Vector)
 * @param width Dimension of stochastic subspace, N_S (number of columns of Vector)
 *
 * @remark TDatabase::ParamDB->COVARIANCE_MATRIX_DO is defined in Database.C
 * @remark TDatabase::ParamDB->COVARIANCE_MATRIX_DO is stored in Row Major form
 */
void CalcCovarianceMatx(double *Vector)
{

	int height = TDatabase::ParamDB->REALISATIONS;
	int width = TDatabase::ParamDB->N_Subspace_Dim;
	// TDatabase::ParamDB->COVARIANCE_MATRIX_DO = new double[width * width]();
	double *phi = new double[width * height](); //Col to Row Major
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			phi[width * i + j] = Vector[height * j + i];
		}
	}

	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, width, width, height, (1.0 / (height - 1)), phi, height, phi, width, 0.0, TDatabase::ParamDB->COVARIANCE_MATRIX_DO, width);


	
}

/**
 * @brief
 *
 * @param Vector Coefficient Matrix
 * @param height Number of Realizations, N_R (number of rows of Vector)
 * @param width Dimension of stochastic subspace, N_S (number of columns of Vector)
 *
 * @remark TDatabase::ParamDB->COVARIANCE_MATRIX_DO is defined in Database.C
 * @remark TDatabase::ParamDB->COVARIANCE_MATRIX_DO is stored in Row Major form
 */
void CalcCoskewnessMatx(double *Vector)
{

	int height = TDatabase::ParamDB->REALISATIONS;
	int width = TDatabase::ParamDB->N_Subspace_Dim;
	double *phi = new double[width * height](); //Col to Row Major
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			phi[width * i + j] = Vector[height * j + i];
		}
	}

	// TDatabase::ParamDB->COSKEWNESS_MATRIX_DO = new double[width * width * width]();
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < width; j++)
		{
			for (int k = 0; k < width; k++)
			{
				for (int p = 0; p < height; p++)
				{

					TDatabase::ParamDB->COSKEWNESS_MATRIX_DO[width * width * k + width * i + j] += ((phi[width * p + i] * phi[width * p + j] * phi[width * p + k]) / (height - 1));
				}
			}
		}
	}
}

lapack_int matInv(double *A, unsigned n)
{
	int ipiv[n + 1];
	lapack_int ret;

	ret = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,
						 n,
						 n,
						 A,
						 n,
						 ipiv);

	if (ret != 0)
		// cout<<"Failure with dgetrf" << endl;
		return ret;

	ret = LAPACKE_dgetri(LAPACK_ROW_MAJOR,
						 n,
						 A,
						 n,
						 ipiv);

	// cout << "Failure with dgetri" << endl;
	return ret;
}

void InvertCov()
{
	int N = TDatabase::ParamDB->N_Subspace_Dim;
	// TDatabase::ParamDB->COVARIANCE_INVERSE_DO = new double[N * N]();

	memcpy(TDatabase::ParamDB->COVARIANCE_INVERSE_DO, TDatabase::ParamDB->COVARIANCE_MATRIX_DO, N * N * SizeOfDouble);
	matInv(TDatabase::ParamDB->COVARIANCE_INVERSE_DO, N);
}

void DO_Mean_RHS(TFESpace2D *Fespace, TFEVectFunct2D *FeVector_Mode,int N_S, double *GlobalRhs_mean,int N_U)
{

	int N_Cells = Fespace->GetN_Cells();
	TCollection *coll = Fespace->GetCollection();
	double* U_Mode = FeVector_Mode->GetValues();
	int lenMode = FeVector_Mode->GetLength();
	// cout << "Length Mode in DO Mean RHS: " << lenMode << endl;
	int lenMean = N_U;


	

	// Get the Global DOF arrays INdex from the FE Space.
	int *GlobalNumbers = Fespace->GetGlobalNumbers();
	int *BeginIndex = Fespace->GetBeginIndex();

	// --- Quadrature Formula Arrays  ------------------//
	int N_Points2;
	double *Weights2, *t1, *t2; // Weights - Quadrature Weights array ; t1  -- Quadrature point ( xi ) in ref coordinate ; t2 -  Quadrature Point ( eta ) in Ref Coordinate
	bool Needs2ndDer[1];
	Needs2ndDer[0] = TRUE;
	double AbsDetjk[MaxN_PointsForNodal2D];
	double X[MaxN_PointsForNodal2D];
	double Y[MaxN_PointsForNodal2D];

	// FE Values Arrays
	double **origvaluesD00; // Shape function values at quadrature Points
	double **origvaluesD10; // Shape Function Derivatives ( x ) at Quadrature Points
	double **origvaluesD01; // Shape Function Derivatives ( y ) at Quadrature Points
	double **origvaluesD20; // Shape Function 2nd Derivatives ( x ) at Quadrature Points
	double **origvaluesD02; // Shape Function 2nd Derivatives ( y ) at Quadrature Points

	for (int cellId = 0; cellId < N_Cells; cellId++)
	{//Cell Loop
		TBaseCell *currentCell = coll->GetCell(cellId);
		// Get the "ID" of Finite Element for the given 2D Element ( Conforming/NonConforming-Order Finite Element : eg : it could be Conforming-2nd order Finite Element )
		FE2D elementId = Fespace->GetFE2D(cellId, currentCell);
		// Get the Class object for that 2d FEM Element , which has all the details like Shape functions , Nodal point locations for that location, Reference Transformation ( Affine , Bilinear )
		TFE2D *element = TFEDatabase2D::GetFE2D(elementId);
		TFEDesc2D *fedesc = element->GetFEDesc2D();
		// Class for basis functions in 2D ( Number of basis functions ), basis function values and Derivatives
		TBaseFunct2D *bf = element->GetBaseFunct2D();
		// Get the Reference Elemet
		BF2DRefElements RefElement = TFEDatabase2D::GetRefElementFromFE2D(elementId);
		// Get the reference Transformation -- Affine Mapping / Bilnea Mapping of Triangle or Quadrilateral
		RefTrans2D referenceTransformation = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(elementId);
		// Get the number of basis functions in the Current Cell ( Number of Local DOF)
		int N_BaseFunct = element->GetN_DOF();
		// Type of Basis Function in 2D
		BaseFunct2D BaseFunct_ID = element->GetBaseFunct2D_ID();

		// get cell measure
		double hK = currentCell->GetDiameter();

		switch (referenceTransformation)
		{
		case QuadBilinear:
		{
			int l = bf->GetPolynomialDegree();									   // Get the Polynomial Degreee  of the basis functions
			QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);		   // Get te ID of Quadrature Formula
			TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
			QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);		   // get the Quadrature points , Weights

			// Set the values on the Reference Cell
			TRefTrans2D *F_K = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
			TFEDatabase2D::SetCellForRefTrans(currentCell, QuadBilinear); // Set the Cell for Current reference Transformation

			// Get Original Coordinates from reference Coordinates and the Determinant of jacobian
			TFEDatabase2D::GetOrigFromRef(QuadBilinear, N_Points2, t1, t2, X, Y, AbsDetjk); // Get the Original Co-orinates for the cell from xi values

			// Get all the original Values from the Referece cell values.
			TFEDatabase2D::GetOrigValues(QuadBilinear, 1, &BaseFunct_ID, N_Points2, t1, t2, QF2, Needs2ndDer);

			// The below are 2D arrays in the form
			// Values[QuadraturePointLocation][ShapeFunction]  i.e, the Value of Shapefunction at all quadrature points for each shape functions
			origvaluesD00 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D00); // Shape Function Values at Quadrature Points
			origvaluesD10 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D10); // Shape Function Derivative Values at Quadrature Points
			origvaluesD01 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D01); // Shape Function Derivative Values at Quadrature Point
			origvaluesD20 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D20); // Shape Function 2nd Derivative Values at Quadrature Point
			origvaluesD02 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D02); // Shape Function 2nd Derivative Values at Quadrature Point
			break;
		}

		case QuadAffin:
		{
			int l = bf->GetPolynomialDegree();									   // Get the Polynomial Degreee  of the basis functions
			QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);		   // Get te ID of Quadrature Formula
			TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
			QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);		   // get the Quadrature points , Weights

			// Set the values on the Reference Cell
			TRefTrans2D *F_K = TFEDatabase2D::GetRefTrans2D(QuadAffin);
			TFEDatabase2D::SetCellForRefTrans(currentCell, QuadAffin); // Set the Cell for Current reference Transformation

			// Get Original Coordinates from reference Coordinates and the Determinant of jacobian
			TFEDatabase2D::GetOrigFromRef(QuadAffin, N_Points2, t1, t2, X, Y, AbsDetjk); // Get the Original Co-orinates for the cell from xi values

			// Get all the original Values from the Referece cell values.
			TFEDatabase2D::GetOrigValues(QuadAffin, 1, &BaseFunct_ID, N_Points2, t1, t2, QF2, Needs2ndDer);

			// The below are 2D arrays in the form
			// Values[QuadraturePointLocation][ShapeFunction]  i.e, the Value of Shapefunction at all quadrature points for each shape functions
			origvaluesD00 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D00); // Shape Function Values at Quadrature Points
			origvaluesD10 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D10); // Shape Function Derivative Values at Quadrature Points
			origvaluesD01 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D01); // Shape Function Derivative Values at Quadrature Point
			origvaluesD20 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D20); // Shape Function 2nd Derivative Values at Quadrature Points
			origvaluesD02 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D02); // Shape Function 2nd Derivative Values at Quadrature Point

			break;
		}

		default:
		{
			cout << " [ERROR] - Error in File : CoeffEqn_DO.C " << endl;
			cout << " Unknown Reftype " << endl;
			cout << " REF TYPE : " << referenceTransformation << endl;
			exit(0);
			break;
		}
		}

		int *DOF = GlobalNumbers + BeginIndex[cellId];
		double val1 = 0;
		double val2 = 0;

		double rhs1[N_BaseFunct];
		double rhs2[N_BaseFunct];
		for (int j = 0; j < N_BaseFunct; j++)
		{
			rhs1[j] = 0;
			rhs2[j] = 0;
		}

		// Get Coefficients b1 and b2
		double *Param[MaxN_QuadPoints_2D];
		double **Coeffs = new double *[MaxN_QuadPoints_2D];
		for (int i = 0; i < MaxN_QuadPoints_2D; i++)
		{
			Coeffs[i] = new double[10]();
		}
		DO_Mean_Equation_Coefficients(N_Points2, X, Y, Param, Coeffs);
	

		for (int a = 0; a < N_S; a++)
		{														   //"a" loop
			double* Mode_Comp1_a = U_Mode + (a * 2 * lenMode);	   // col Major
			double* Mode_Comp2_a = U_Mode + ((a * 2 + 1) * lenMode); // col Major
																   //  double* phi_Array_a = Phi_Array + a*lenMode;??

			double U1_Mode_a[N_Points2];
			double U1x_Mode_a[N_Points2];
			double U1y_Mode_a[N_Points2];

			double U2_Mode_a[N_Points2];
			double U2x_Mode_a[N_Points2];
			double U2y_Mode_a[N_Points2];

			for (int quadPt = 0; quadPt < N_Points2; quadPt++) // Initialize
			{
				// C_i[quadPt] = 0;
				U1_Mode_a[quadPt] = 0;
				U1x_Mode_a[quadPt] = 0;
				U1y_Mode_a[quadPt] = 0;

				U2_Mode_a[quadPt] = 0;
				U2x_Mode_a[quadPt] = 0;
				U2y_Mode_a[quadPt] = 0;
			}
			

			// Obtain all values for C_a
			for (int quadPt = 0; quadPt < N_Points2; quadPt++)
			{
				for (int j = 0; j < N_BaseFunct; j++)
				{
					int globDOF = DOF[j];
					U1_Mode_a[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp1_a[globDOF];
					U1x_Mode_a[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp1_a[globDOF];
					U1y_Mode_a[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp1_a[globDOF];

					U2_Mode_a[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp2_a[globDOF];
					U2x_Mode_a[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp2_a[globDOF];
					U2y_Mode_a[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp2_a[globDOF];
				}
			}
			for (int b = 0; b < N_S; b++)
			{ //"b" loop
				val1 = 0;
				val2 = 0;
				double* Mode_Comp1_b = U_Mode + (b * 2 * lenMode);	   // col Major
				double* Mode_Comp2_b = U_Mode + ((b * 2 + 1) * lenMode); // col Major

				double U1_Mode_b[N_Points2];
				double U1x_Mode_b[N_Points2];
				double U1y_Mode_b[N_Points2];

				double U2_Mode_b[N_Points2];
				double U2x_Mode_b[N_Points2];
				double U2y_Mode_b[N_Points2];

				for (int quadPt = 0; quadPt < N_Points2; quadPt++) // Initialize
				{
					// C_i[quadPt] = 0;
					U1_Mode_b[quadPt] = 0;
					U1x_Mode_b[quadPt] = 0;
					U1y_Mode_b[quadPt] = 0;

					U2_Mode_b[quadPt] = 0;
					U2x_Mode_b[quadPt] = 0;
					U2y_Mode_b[quadPt] = 0;
				}

				for (int quadPt = 0; quadPt < N_Points2; quadPt++)
				{
					for (int j = 0; j < N_BaseFunct; j++)
					{
						int globDOF = DOF[j];
						U1_Mode_b[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp1_b[globDOF];
						U1x_Mode_b[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp1_b[globDOF];
						U1y_Mode_b[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp1_b[globDOF];

						U2_Mode_b[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp2_b[globDOF];
						U2x_Mode_b[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp2_b[globDOF];
						U2y_Mode_b[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp2_b[globDOF];
					}
				}

				for (int qdpt = 0; qdpt < N_Points2; qdpt++)
				{
					double Mult = Weights2[qdpt] * AbsDetjk[qdpt];
					double *orgD00 = origvaluesD00[qdpt];

					val1 += -1.0 * (U1_Mode_a[qdpt] * U1x_Mode_b[qdpt] + U2_Mode_a[qdpt] * U1y_Mode_b[qdpt])* TDatabase::ParamDB->COVARIANCE_MATRIX_DO[a * N_S + b] * Mult;

					val2 += -1.0 * (U1_Mode_a[qdpt] * U2x_Mode_b[qdpt] + U2_Mode_a[qdpt] * U2y_Mode_b[qdpt])* TDatabase::ParamDB->COVARIANCE_MATRIX_DO[a * N_S + b] * Mult;

					for (int j = 0; j < N_BaseFunct; j++)
					{
						rhs1[j] += val1 * orgD00[j];
						rhs2[j] += val2 * orgD00[j]; // * Mult;
					}
				}

			} //"b" loop

		} //"a" loop


		for (int j = 0; j < N_BaseFunct; j++)
		{
			int GlobalDOF = DOF[j];
			GlobalRhs_mean[GlobalDOF] += rhs1[j];
			GlobalRhs_mean[GlobalDOF+lenMean] += rhs2[j];
		}
		// --
	}//cell loop

	
}

//======================================================================
// ************************** Mode RHS********************************//
//======================================================================
void DO_Mode_RHS(TFESpace2D *Fespace, TFEVectFunct2D *FeVector_Mean, TFEVectFunct2D *FeVector_Mode, int N_S, double *GlobalRhs_mode, int i_index)
{

	int N_Cells = Fespace->GetN_Cells();
	TCollection *coll = Fespace->GetCollection();

	// Get the Global DOF arrays INdex from the FE Space.
	int *GlobalNumbers = Fespace->GetGlobalNumbers();
	int *BeginIndex = Fespace->GetBeginIndex();

	// --- Quadrature Formula Arrays  ------------------//
	int N_Points2;
	double *Weights2, *t1, *t2; // Weights - Quadrature Weights array ; t1  -- Quadrature point ( xi ) in ref coordinate ; t2 -  Quadrature Point ( eta ) in Ref Coordinate
	bool Needs2ndDer[1];
	Needs2ndDer[0] = TRUE;
	double AbsDetjk[MaxN_PointsForNodal2D];
	double X[MaxN_PointsForNodal2D];
	double Y[MaxN_PointsForNodal2D];

	// FE Values Arrays
	double **origvaluesD00; // Shape function values at quadrature Points
	double **origvaluesD10; // Shape Function Derivatives ( x ) at Quadrature Points
	double **origvaluesD01; // Shape Function Derivatives ( y ) at Quadrature Points
	double **origvaluesD20; // Shape Function 2nd Derivatives ( x ) at Quadrature Points
	double **origvaluesD02; // Shape Function 2nd Derivatives ( y ) at Quadrature Points

	double *U_Mode = FeVector_Mode->GetValues();
	int lenMode = FeVector_Mode->GetLength();

	double *U_Mean = FeVector_Mean->GetValues();
	int lenMean = FeVector_Mean->GetLength();

	double *Mean_Comp1 = U_Mean;
	double *Mean_Comp2 = U_Mean + lenMean; //

	double *Mode_Comp1_i = U_Mode + i_index * 2 * lenMode;		 // col Major
	double *Mode_Comp2_i = U_Mode + (i_index * 2 + 1) * lenMode; // col Major

	double val1 = 0;
	double val2 = 0;

	for (int cellId = 0; cellId < N_Cells; cellId++)
	{
		TBaseCell *currentCell = coll->GetCell(cellId);
		// Get the "ID" of Finite Element for the given 2D Element ( Conforming/NonConforming-Order Finite Element : eg : it could be Conforming-2nd order Finite Element )
		FE2D elementId = Fespace->GetFE2D(cellId, currentCell);
		// Get the Class object for that 2d FEM Element , which has all the details like Shape functions , Nodal point locations for that location, Reference Transformation ( Affine , Bilinear )
		TFE2D *element = TFEDatabase2D::GetFE2D(elementId);
		TFEDesc2D *fedesc = element->GetFEDesc2D();
		// Class for basis functions in 2D ( Number of basis functions ), basis function values and Derivatives
		TBaseFunct2D *bf = element->GetBaseFunct2D();
		// Get the Reference Elemet
		BF2DRefElements RefElement = TFEDatabase2D::GetRefElementFromFE2D(elementId);
		// Get the reference Transformation -- Affine Mapping / Bilnea Mapping of Triangle or Quadrilateral
		RefTrans2D referenceTransformation = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(elementId);
		// Get the number of basis functions in the Current Cell ( Number of Local DOF)
		int N_BaseFunct = element->GetN_DOF();
		// Type of Basis Function in 2D
		BaseFunct2D BaseFunct_ID = element->GetBaseFunct2D_ID();

		// get cell measure
		double hK = currentCell->GetDiameter();

		switch (referenceTransformation)
		{
		case QuadBilinear:
		{
			int l = bf->GetPolynomialDegree();									   // Get the Polynomial Degreee  of the basis functions
			QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);		   // Get te ID of Quadrature Formula
			TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
			QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);		   // get the Quadrature points , Weights

			// Set the values on the Reference Cell
			TRefTrans2D *F_K = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
			TFEDatabase2D::SetCellForRefTrans(currentCell, QuadBilinear); // Set the Cell for Current reference Transformation

			// Get Original Coordinates from reference Coordinates and the Determinant of jacobian
			TFEDatabase2D::GetOrigFromRef(QuadBilinear, N_Points2, t1, t2, X, Y, AbsDetjk); // Get the Original Co-orinates for the cell from xi values

			// Get all the original Values from the Referece cell values.
			TFEDatabase2D::GetOrigValues(QuadBilinear, 1, &BaseFunct_ID, N_Points2, t1, t2, QF2, Needs2ndDer);

			// The below are 2D arrays in the form
			// Values[QuadraturePointLocation][ShapeFunction]  i.e, the Value of Shapefunction at all quadrature points for each shape functions
			origvaluesD00 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D00); // Shape Function Values at Quadrature Points
			origvaluesD10 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D10); // Shape Function Derivative Values at Quadrature Points
			origvaluesD01 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D01); // Shape Function Derivative Values at Quadrature Point
			origvaluesD20 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D20); // Shape Function 2nd Derivative Values at Quadrature Point
			origvaluesD02 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D02); // Shape Function 2nd Derivative Values at Quadrature Point
			break;
		}

		case QuadAffin:
		{
			int l = bf->GetPolynomialDegree();									   // Get the Polynomial Degreee  of the basis functions
			QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);		   // Get te ID of Quadrature Formula
			TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
			QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);		   // get the Quadrature points , Weights

			// Set the values on the Reference Cell
			TRefTrans2D *F_K = TFEDatabase2D::GetRefTrans2D(QuadAffin);
			TFEDatabase2D::SetCellForRefTrans(currentCell, QuadAffin); // Set the Cell for Current reference Transformation

			// Get Original Coordinates from reference Coordinates and the Determinant of jacobian
			TFEDatabase2D::GetOrigFromRef(QuadAffin, N_Points2, t1, t2, X, Y, AbsDetjk); // Get the Original Co-orinates for the cell from xi values

			// Get all the original Values from the Referece cell values.
			TFEDatabase2D::GetOrigValues(QuadAffin, 1, &BaseFunct_ID, N_Points2, t1, t2, QF2, Needs2ndDer);

			// The below are 2D arrays in the form
			// Values[QuadraturePointLocation][ShapeFunction]  i.e, the Value of Shapefunction at all quadrature points for each shape functions
			origvaluesD00 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D00); // Shape Function Values at Quadrature Points
			origvaluesD10 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D10); // Shape Function Derivative Values at Quadrature Points
			origvaluesD01 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D01); // Shape Function Derivative Values at Quadrature Point
			origvaluesD20 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D20); // Shape Function 2nd Derivative Values at Quadrature Points
			origvaluesD02 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D02); // Shape Function 2nd Derivative Values at Quadrature Point

			break;
		}

		default:
		{
			cout << " [ERROR] - Error in File : CoeffEqn_DO.C " << endl;
			cout << " Unknown Reftype " << endl;
			cout << " REF TYPE : " << referenceTransformation << endl;
			exit(0);
			break;
		}
		}

		int *DOF = GlobalNumbers + BeginIndex[cellId];
		double val = 0;

		// Get Coefficients b1 and b2
		double *Param[MaxN_QuadPoints_2D];
		double **Coeffs = new double *[MaxN_QuadPoints_2D];
		for (int i = 0; i < MaxN_QuadPoints_2D; i++)
		{
			Coeffs[i] = new double[10]();
		}

		double U1_Mode_i[N_Points2];
		double U1x_Mode_i[N_Points2];
		double U1y_Mode_i[N_Points2];
		double U1xx_Mode_i[N_Points2];
		double U1yy_Mode_i[N_Points2];

		double U2_Mode_i[N_Points2];
		double U2x_Mode_i[N_Points2];
		double U2y_Mode_i[N_Points2];
		double U2xx_Mode_i[N_Points2];
		double U2yy_Mode_i[N_Points2];

		double U1_Mean[N_Points2];
		double U1x_Mean[N_Points2];
		double U1y_Mean[N_Points2];

		double U2_Mean[N_Points2];
		double U2x_Mean[N_Points2];
		double U2y_Mean[N_Points2];

		for (int quadPt = 0; quadPt < N_Points2; quadPt++) // Initialize
		{
			U1_Mode_i[quadPt] = 0;
			U1x_Mode_i[quadPt] = 0;
			U1y_Mode_i[quadPt] = 0;
			U1xx_Mode_i[quadPt] = 0;
			U1yy_Mode_i[quadPt] = 0;

			U2_Mode_i[quadPt] = 0;
			U2x_Mode_i[quadPt] = 0;
			U2y_Mode_i[quadPt] = 0;
			U2xx_Mode_i[quadPt] = 0;
			U2yy_Mode_i[quadPt] = 0;

			U1_Mean[quadPt] = 0;
			U1x_Mean[quadPt] = 0;
			U1y_Mean[quadPt] = 0;

			U2_Mean[quadPt] = 0;
			U2x_Mean[quadPt] = 0;
			U2y_Mean[quadPt] = 0;
		}

		for (int quadPt = 0; quadPt < N_Points2; quadPt++)
		{
			for (int j = 0; j < N_BaseFunct; j++)
			{
				int globDOF = DOF[j];
				U1_Mode_i[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp1_i[globDOF];
				U1x_Mode_i[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp1_i[globDOF];
				U1y_Mode_i[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp1_i[globDOF];
				U1xx_Mode_i[quadPt] += origvaluesD20[quadPt][j] * Mode_Comp1_i[globDOF];
				U1yy_Mode_i[quadPt] += origvaluesD02[quadPt][j] * Mode_Comp1_i[globDOF];

				U2_Mode_i[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp2_i[globDOF];
				U2x_Mode_i[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp2_i[globDOF];
				U2y_Mode_i[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp2_i[globDOF];
				U2xx_Mode_i[quadPt] += origvaluesD20[quadPt][j] * Mode_Comp2_i[globDOF];
				U2yy_Mode_i[quadPt] += origvaluesD02[quadPt][j] * Mode_Comp2_i[globDOF];

				U1_Mean[quadPt] += origvaluesD00[quadPt][j] * Mean_Comp1[globDOF];
				U1x_Mean[quadPt] += origvaluesD10[quadPt][j] * Mean_Comp1[globDOF];
				U1y_Mean[quadPt] += origvaluesD01[quadPt][j] * Mean_Comp1[globDOF];

				U2_Mean[quadPt] += origvaluesD00[quadPt][j] * Mean_Comp2[globDOF];
				U2x_Mean[quadPt] += origvaluesD10[quadPt][j] * Mean_Comp2[globDOF];
				U2y_Mean[quadPt] += origvaluesD01[quadPt][j] * Mean_Comp2[globDOF];
			}
		}

		// Save Values of C at all quadrature points for I component

		double rhs1[N_BaseFunct];
		double rhs2[N_BaseFunct];
		for (int j = 0; j < N_BaseFunct; j++)
		{
			rhs1[j] = 0;
			rhs2[j] = 0;
		}

		for (int a = 0; a < N_S; a++)
		{ //"a" loop

			double *Mode_Comp1_a = U_Mode + a * 2 * lenMode;	   // col Major
			double *Mode_Comp2_a = U_Mode + (a * 2 + 1) * lenMode; // col Major

			double U1_Mode_a[N_Points2];
			double U1x_Mode_a[N_Points2];
			double U1y_Mode_a[N_Points2];
			double U1xx_Mode_a[N_Points2];
			double U1yy_Mode_a[N_Points2];

			double U2_Mode_a[N_Points2];
			double U2x_Mode_a[N_Points2];
			double U2y_Mode_a[N_Points2];
			double U2xx_Mode_a[N_Points2];
			double U2yy_Mode_a[N_Points2];

			for (int quadPt = 0; quadPt < N_Points2; quadPt++) // Initialize
			{
				// C_i[quadPt] = 0;
				U1_Mode_a[quadPt] = 0;
				U1x_Mode_a[quadPt] = 0;
				U1y_Mode_a[quadPt] = 0;

				U2_Mode_a[quadPt] = 0;
				U2x_Mode_a[quadPt] = 0;
				U2y_Mode_a[quadPt] = 0;
			}
			// for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++) C_x_a[quadPt] = 0;
			// for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++) C_y_a[quadPt] = 0;

			// Obtain all values for U_a
			for (int quadPt = 0; quadPt < N_Points2; quadPt++)
			{
				for (int j = 0; j < N_BaseFunct; j++)
				{
					int globDOF = DOF[j];
					U1_Mode_a[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp1_a[globDOF];
					U1x_Mode_a[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp1_a[globDOF];
					U1y_Mode_a[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp1_a[globDOF];

					U2_Mode_a[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp2_a[globDOF];
					U2x_Mode_a[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp2_a[globDOF];
					U2y_Mode_a[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp2_a[globDOF];
				}
			}

			for (int b = 0; b < N_S; b++)
			{ //"b" loop starts

				double *Mode_Comp1_b = U_Mode + b * 2 * lenMode;	   // col Major
				double *Mode_Comp2_b = U_Mode + (b * 2 + 1) * lenMode; // col Major

				double U1_Mode_b[N_Points2];
				double U1x_Mode_b[N_Points2];
				double U1y_Mode_b[N_Points2];

				double U2_Mode_b[N_Points2];
				double U2x_Mode_b[N_Points2];
				double U2y_Mode_b[N_Points2];

				for (int quadPt = 0; quadPt < N_Points2; quadPt++) // Initialize
				{
					U1_Mode_b[quadPt] = 0;
					U1x_Mode_b[quadPt] = 0;
					U1y_Mode_b[quadPt] = 0;

					U2_Mode_b[quadPt] = 0;
					U2x_Mode_b[quadPt] = 0;
					U2y_Mode_b[quadPt] = 0;
				}

				for (int quadPt = 0; quadPt < N_Points2; quadPt++)
				{
					for (int j = 0; j < N_BaseFunct; j++)
					{
						int globDOF = DOF[j];
						U1_Mode_b[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp1_b[globDOF];
						U1x_Mode_b[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp1_b[globDOF];
						U1y_Mode_b[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp1_b[globDOF];

						U2_Mode_b[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp2_b[globDOF];
						U2x_Mode_b[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp2_b[globDOF];
						U2y_Mode_b[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp2_b[globDOF];
					}
				}

				for (int p = 0; p < N_S; p++)
				{ //"p" loop starts
					val1 = 0;
					val2 = 0;
					double *Mode_Comp1_p = U_Mode + p * 2 * lenMode;	   // col Major
					double *Mode_Comp2_p = U_Mode + (p * 2 + 1) * lenMode; // col Major

					double U1_Mode_p[N_Points2];
					double U1x_Mode_p[N_Points2];
					double U1y_Mode_p[N_Points2];

					double U2_Mode_p[N_Points2];
					double U2x_Mode_p[N_Points2];
					double U2y_Mode_p[N_Points2];

					for (int quadPt = 0; quadPt < N_Points2; quadPt++) // Initialize
					{
						U1_Mode_p[quadPt] = 0;
						U1x_Mode_p[quadPt] = 0;
						U1y_Mode_p[quadPt] = 0;

						U2_Mode_p[quadPt] = 0;
						U2x_Mode_p[quadPt] = 0;
						U2y_Mode_p[quadPt] = 0;
					}

					for (int quadPt = 0; quadPt < N_Points2; quadPt++)
					{
						for (int j = 0; j < N_BaseFunct; j++)
						{
							int globDOF = DOF[j];
							U1_Mode_p[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp1_p[globDOF];
							U1x_Mode_p[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp1_p[globDOF];
							U1y_Mode_p[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp1_p[globDOF];

							U2_Mode_p[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp2_p[globDOF];
							U2x_Mode_p[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp2_p[globDOF];
							U2y_Mode_p[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp2_p[globDOF];
						}
					}

					// INner Quadrature Loop
					for (int quadPt = 0; quadPt < N_Points2; quadPt++)
					{
						double Mult = Weights2[quadPt] * AbsDetjk[quadPt];
						double *orgD00 = origvaluesD00[quadPt];
						double *orgD10 = origvaluesD10[quadPt];
						double *orgD01 = origvaluesD01[quadPt];
						double *orgD20 = origvaluesD20[quadPt];
						double *orgD02 = origvaluesD02[quadPt];

						double nu = Coeffs[quadPt][0]; // nu

						for (int c = 0; c < N_S; c++)
						{ //"c"loop
							for (int qdpt = 0; qdpt < N_Points2; qdpt++)
							{
								double Mult = Weights2[qdpt] * AbsDetjk[qdpt];

								val1 += -1.0 * (TDatabase::ParamDB->COVARIANCE_INVERSE_DO[N_S * c + i_index] * TDatabase::ParamDB->COSKEWNESS_MATRIX_DO[N_S * N_S * b + N_S * c + a]) * (U1_Mode_a[qdpt] * U1x_Mode_b[qdpt] + U2_Mode_a[qdpt] * U1y_Mode_b[qdpt]) * U1_Mode_p[qdpt] * Mult;

								val2 += -1.0 * (TDatabase::ParamDB->COVARIANCE_INVERSE_DO[N_S * c + i_index] * TDatabase::ParamDB->COSKEWNESS_MATRIX_DO[N_S * N_S * b + N_S * c + a]) * (U1_Mode_a[qdpt] * U2x_Mode_b[qdpt] + U2_Mode_a[qdpt] * U2y_Mode_b[qdpt]) * U2_Mode_p[qdpt] * Mult;
							}

						} //"c" loop ends

						for (int qdpt = 0; qdpt < N_Points2; qdpt++)
						{
							val1 += (-1.0 * (U1_Mode_i[qdpt] * U1x_Mean[qdpt] + U2_Mode_a[qdpt] * U1y_Mean[qdpt])) * U1_Mode_p[qdpt] * Mult;

							val1 += (-1.0 * (U1_Mean[qdpt] * U1x_Mode_i[qdpt] + U2_Mean[qdpt] * U1y_Mode_i[qdpt])) * U1_Mode_p[qdpt] * Mult;

							val1 += nu * (U1xx_Mode_i[qdpt] + U1yy_Mode_i[qdpt]) * U1_Mode_p[qdpt] * Mult;

							val2 += (-1.0 * (U1_Mode_i[qdpt] * U2x_Mean[qdpt] + U2_Mode_a[qdpt] * U2y_Mean[qdpt])) * U2_Mode_p[qdpt] * Mult;

							val2 += (-1.0 * (U1_Mean[qdpt] * U2x_Mode_i[qdpt] + U2_Mean[qdpt] * U2y_Mode_i[qdpt])) * U2_Mode_p[qdpt] * Mult;

							val2 += nu * (U2xx_Mode_i[qdpt] + U2yy_Mode_i[qdpt]) * U2_Mode_p[qdpt] * Mult;
						}

						val1 *= -1.0 * U1_Mode_p[quadPt]; // This is Final "f"
						val2 *= -1.0 * U2_Mode_p[quadPt]; // This is Final "f"

						for (int j = 0; j < N_BaseFunct; j++)
						{
							rhs1[j] += val1 * orgD00[j]; // * Mult;
							rhs2[j] += val2 * orgD00[j]; // * Mult;
						}

					} // Inner Quadrature Loop for p

				} //"p" loop ends
				val1 = 0;
				val2 = 0;
				for (int c = 0; c < N_S; c++)
				{ //"c" loop
					for (int qdpt = 0; qdpt < N_Points2; qdpt++)
					{
						double Mult = Weights2[qdpt] * AbsDetjk[qdpt];
						double *orgD00 = origvaluesD00[qdpt];
						double *orgD10 = origvaluesD10[qdpt];
						double *orgD01 = origvaluesD01[qdpt];
						double *orgD20 = origvaluesD20[qdpt];
						double *orgD02 = origvaluesD02[qdpt];

						val1 += -1.0 * (TDatabase::ParamDB->COVARIANCE_INVERSE_DO[N_S * c + i_index] * TDatabase::ParamDB->COSKEWNESS_MATRIX_DO[N_S * N_S * b + N_S * c + a]) * (U1_Mode_a[qdpt] * U1x_Mode_b[qdpt] + U2_Mode_a[qdpt] * U1y_Mode_b[qdpt]) * Mult;

						val2 += -1.0 * (TDatabase::ParamDB->COVARIANCE_INVERSE_DO[N_S * c + i_index] * TDatabase::ParamDB->COSKEWNESS_MATRIX_DO[N_S * N_S * b + N_S * c + a]) * (U1_Mode_a[qdpt] * U2x_Mode_b[qdpt] + U2_Mode_a[qdpt] * U2y_Mode_b[qdpt]) * Mult;

						for (int j = 0; j < N_BaseFunct; j++)
						{
							rhs1[j] += val1 * orgD00[j]; // * Mult;
							rhs2[j] += val2 * orgD00[j]; // * Mult;
						}

					} // Internal quadrature point loop
				}	  //"c" loop ends here

			} //"b" loop ends

		} //"a" loop ends
		val1 = 0;
		val2 = 0;
		for (int qdpt = 0; qdpt < N_Points2; qdpt++)
		{
			double Mult = Weights2[qdpt] * AbsDetjk[qdpt];
			double *orgD00 = origvaluesD00[qdpt];
			double *orgD10 = origvaluesD10[qdpt];
			double *orgD01 = origvaluesD01[qdpt];
			double *orgD20 = origvaluesD20[qdpt];
			double *orgD02 = origvaluesD02[qdpt];

			val1 += -1.0 * (U1_Mode_i[qdpt] * U1x_Mean[qdpt] + U2_Mode_i[qdpt] * U1y_Mean[qdpt]) * Mult;
			val1 += -1.0 * (U1_Mean[qdpt] * U1x_Mode_i[qdpt] + U2_Mode_i[qdpt] * U1y_Mode_i[qdpt]) * Mult;
			val1 += (U1_Mode_i[qdpt] * U1x_Mode_i[qdpt] + U2_Mode_i[qdpt] * U1y_Mode_i[qdpt]) * Mult;

			val2 += -1.0 * (U1_Mode_i[qdpt] * U2x_Mean[qdpt] + U2_Mode_i[qdpt] * U2y_Mean[qdpt]) * Mult;
			val2 += -1.0 * (U1_Mean[qdpt] * U2x_Mode_i[qdpt] + U2_Mode_i[qdpt] * U2y_Mode_i[qdpt]) * Mult;
			val2 += (U1_Mode_i[qdpt] * U2x_Mode_i[qdpt] + U2_Mode_i[qdpt] * U2y_Mode_i[qdpt]) * Mult;

			for (int j = 0; j < N_BaseFunct; j++)
			{
				rhs1[j] += val1 * orgD00[j]; // * Mult;
				rhs2[j] += val2 * orgD00[j]; // * Mult;
			}

		} // Internal quadrature point loop

		for (int j = 0; j < N_BaseFunct; j++)
		{
			int GlobalDOF = DOF[j];
			GlobalRhs_mode[GlobalDOF] += rhs1[j];
			GlobalRhs_mode[GlobalDOF + lenMode] += rhs2[j];
		}
		// --
	}
}

//======================================================================
// ************************** Coefficient ********************************//
//======================================================================
/**
 * @brief
 *
 * @param Fespace
 * @param FeVector_Mode
 * @param FEVector_Phi
 * @param FeVector_Mean
 * @param N_S
 * @param i_index
 * @param N_R
 */
void DO_CoEfficient(TFESpace2D *Fespace, TFEVectFunct2D *FeVector_Mode, TFEVectFunct2D *FEVector_Phi, TFEVectFunct2D *FeVector_Mean, int N_S, int i_index, int N_R)
{

	int N_Cells = Fespace->GetN_Cells();
	TCollection *coll = Fespace->GetCollection();

	// Get the Global DOF arrays INdex from the FE Space.
	int *GlobalNumbers = Fespace->GetGlobalNumbers();
	int *BeginIndex = Fespace->GetBeginIndex();

	// --- Quadrature Formula Arrays  ------------------//
	int N_Points2;
	double *Weights2, *t1, *t2; // Weights - Quadrature Weights array ; t1  -- Quadrature point ( xi ) in ref coordinate ; t2 -  Quadrature Point ( eta ) in Ref Coordinate
	bool Needs2ndDer[1];
	Needs2ndDer[0] = TRUE;
	double AbsDetjk[MaxN_PointsForNodal2D];
	double X[MaxN_PointsForNodal2D];
	double Y[MaxN_PointsForNodal2D];

	// FE Values Arrays
	double **origvaluesD00; // Shape function values at quadrature Points
	double **origvaluesD10; // Shape Function Derivatives ( x ) at Quadrature Points
	double **origvaluesD01; // Shape Function Derivatives ( y ) at Quadrature Points
	double **origvaluesD20; // Shape Function 2nd Derivatives ( x ) at Quadrature Points
	double **origvaluesD02; // Shape Function 2nd Derivatives ( y ) at Quadrature Points

	double *U_Mode = FeVector_Mode->GetValues();
	int lenMode = FeVector_Mode->GetLength();

	// Check
	double *U_Mean = FeVector_Mean->GetValues();
	int lenMean = FeVector_Mean->GetLength();
	double *Mean_Comp1 = U_Mean;
	double *Mean_Comp2 = U_Mean + lenMean; //

	double *Phi_Array = FEVector_Phi->GetValues();
	int lenPhi = FEVector_Phi->GetLength();

	double *Phi_Old = new double[N_R * N_S]();
	memcpy(Phi_Old, Phi_Array, N_R * N_S * SizeOfDouble);

	double *phi_New = new double[lenPhi]();

	// Check
	double *Mode_Comp1_i = U_Mode + i_index * 2 * lenMode;		 // col Major
	double *Mode_Comp2_i = U_Mode + (i_index * 2 + 1) * lenMode; // col Major
	// 2nd Component??

	double *phi_Array_i = Phi_Array + i_index * lenPhi;
	double *phi_Old_i = Phi_Old + i_index * lenPhi;

	//   double* phi_New = new double[lenPhi]();

	for (int cellId = 0; cellId < N_Cells; cellId++)
	{ // cell loop
		TBaseCell *currentCell = coll->GetCell(cellId);
		// Get the "ID" of Finite Element for the given 2D Element ( Conforming/NonConforming-Order Finite Element : eg : it could be Conforming-2nd order Finite Element )
		FE2D elementId = Fespace->GetFE2D(cellId, currentCell);
		// Get the Class object for that 2d FEM Element , which has all the details like Shape functions , Nodal point locations for that location, Reference Transformation ( Affine , Bilinear )
		TFE2D *element = TFEDatabase2D::GetFE2D(elementId);
		TFEDesc2D *fedesc = element->GetFEDesc2D();
		// Class for basis functions in 2D ( Number of basis functions ), basis function values and Derivatives
		TBaseFunct2D *bf = element->GetBaseFunct2D();
		// Get the Reference Element
		BF2DRefElements RefElement = TFEDatabase2D::GetRefElementFromFE2D(elementId);
		// Get the reference Transformation -- Affine Mapping / Bilnear Mapping of Triangle or Quadrilateral
		RefTrans2D referenceTransformation = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(elementId);
		// Get the number of basis functions in the Current Cell ( Number of Local DOF)
		int N_BaseFunct = element->GetN_DOF();
		// Type of Basis Function in 2D
		BaseFunct2D BaseFunct_ID = element->GetBaseFunct2D_ID();

		// get cell measure
		double hK = currentCell->GetDiameter();

		switch (referenceTransformation)
		{
		case QuadBilinear:
		{
			int l = bf->GetPolynomialDegree();									   // Get the Polynomial Degreee  of the basis functions
			QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);		   // Get te ID of Quadrature Formula
			TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
			QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);		   // get the Quadrature points , Weights

			// Set the values on the Reference Cell
			TRefTrans2D *F_K = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
			TFEDatabase2D::SetCellForRefTrans(currentCell, QuadBilinear); // Set the Cell for Current reference Transformation

			// Get Original Coordinates from reference Coordinates and the Determinant of jacobian
			TFEDatabase2D::GetOrigFromRef(QuadBilinear, N_Points2, t1, t2, X, Y, AbsDetjk); // Get the Original Co-orinates for the cell from xi values

			// Get all the original Values from the Referece cell values.
			TFEDatabase2D::GetOrigValues(QuadBilinear, 1, &BaseFunct_ID, N_Points2, t1, t2, QF2, Needs2ndDer);

			// The below are 2D arrays in the form
			// Values[QuadraturePointLocation][ShapeFunction]  i.e, the Value of Shapefunction at all quadrature points for each shape functions
			origvaluesD00 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D00); // Shape Function Values at Quadrature Points
			origvaluesD10 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D10); // Shape Function Derivative Values at Quadrature Points
			origvaluesD01 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D01); // Shape Function Derivative Values at Quadrature Point
			origvaluesD20 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D20); // Shape Function 2nd Derivative Values at Quadrature Point
			origvaluesD02 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D02); // Shape Function 2nd Derivative Values at Quadrature Point
			break;
		}

		case QuadAffin:
		{
			int l = bf->GetPolynomialDegree();									   // Get the Polynomial Degreee  of the basis functions
			QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);		   // Get te ID of Quadrature Formula
			TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
			QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);		   // get the Quadrature points , Weights

			// Set the values on the Reference Cell
			TRefTrans2D *F_K = TFEDatabase2D::GetRefTrans2D(QuadAffin);
			TFEDatabase2D::SetCellForRefTrans(currentCell, QuadAffin); // Set the Cell for Current reference Transformation

			// Get Original Coordinates from reference Coordinates and the Determinant of jacobian
			TFEDatabase2D::GetOrigFromRef(QuadAffin, N_Points2, t1, t2, X, Y, AbsDetjk); // Get the Original Co-orinates for the cell from xi values

			// Get all the original Values from the Referece cell values.
			TFEDatabase2D::GetOrigValues(QuadAffin, 1, &BaseFunct_ID, N_Points2, t1, t2, QF2, Needs2ndDer);

			// The below are 2D arrays in the form
			// Values[QuadraturePointLocation][ShapeFunction]  i.e, the Value of Shapefunction at all quadrature points for each shape functions
			origvaluesD00 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D00); // Shape Function Values at Quadrature Points
			origvaluesD10 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D10); // Shape Function Derivative Values at Quadrature Points
			origvaluesD01 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D01); // Shape Function Derivative Values at Quadrature Point
			origvaluesD20 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D20); // Shape Function 2nd Derivative Values at Quadrature Points
			origvaluesD02 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D02); // Shape Function 2nd Derivative Values at Quadrature Point

			break;
		}

		default:
		{
			cout << " [ERROR] - Error in File : CoeffEqn_DO.C " << endl;
			cout << " Unknown Reftype " << endl;
			cout << " REF TYPE : " << referenceTransformation << endl;
			exit(0);
			break;
		}
		}

		int *DOF = GlobalNumbers + BeginIndex[cellId];
		double val = 0;

		// Get Coefficients b1 and b2
		double *Param[MaxN_QuadPoints_2D];
		double **Coeffs = new double *[MaxN_QuadPoints_2D];
		for (int i = 0; i < MaxN_QuadPoints_2D; i++)
		{
			Coeffs[i] = new double[10]();
		}

		// Save Values of C at all quadrature points for I component
		double C_i[N_Points2];
		double C_x_i[N_Points2];
		double C_y_i[N_Points2];

		double U1_Mode_i[N_Points2];
		double U1x_Mode_i[N_Points2];
		double U1y_Mode_i[N_Points2];

		double U2_Mode_i[N_Points2];
		double U2x_Mode_i[N_Points2];
		double U2y_Mode_i[N_Points2];

		double U1_Mean[N_Points2];
		double U1x_Mean[N_Points2];
		double U1y_Mean[N_Points2];

		double U2_Mean[N_Points2];
		double U2x_Mean[N_Points2];
		double U2y_Mean[N_Points2];

		for (int quadPt = 0; quadPt < N_Points2; quadPt++) // Initialize
		{
			// C_i[quadPt] = 0;
			U1_Mode_i[quadPt] = 0;
			U1x_Mode_i[quadPt] = 0;
			U1y_Mode_i[quadPt] = 0;

			U2_Mode_i[quadPt] = 0;
			U2x_Mode_i[quadPt] = 0;
			U2y_Mode_i[quadPt] = 0;

			U1_Mean[quadPt] = 0;
			U1x_Mean[quadPt] = 0;
			U1y_Mean[quadPt] = 0;

			U2_Mean[quadPt] = 0;
			U2x_Mean[quadPt] = 0;
			U2y_Mean[quadPt] = 0;
		}
		// for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++) C_x_i[quadPt] = 0;
		// for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++) C_y_i[quadPt] = 0;

		for (int quadPt = 0; quadPt < N_Points2; quadPt++)
		{
			for (int j = 0; j < N_BaseFunct; j++)
			{
				int globDOF = DOF[j];
				U1_Mode_i[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp1_i[globDOF];
				U1x_Mode_i[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp1_i[globDOF];
				U1y_Mode_i[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp1_i[globDOF];

				U2_Mode_i[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp2_i[globDOF];
				U2x_Mode_i[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp2_i[globDOF];
				U2y_Mode_i[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp2_i[globDOF];

				U1_Mean[quadPt] += origvaluesD00[quadPt][j] * Mean_Comp1[globDOF];
				U1x_Mean[quadPt] += origvaluesD10[quadPt][j] * Mean_Comp1[globDOF];
				U1y_Mean[quadPt] += origvaluesD01[quadPt][j] * Mean_Comp1[globDOF];

				U2_Mean[quadPt] += origvaluesD00[quadPt][j] * Mean_Comp2[globDOF];
				U2x_Mean[quadPt] += origvaluesD10[quadPt][j] * Mean_Comp2[globDOF];
				U2y_Mean[quadPt] += origvaluesD01[quadPt][j] * Mean_Comp2[globDOF];

				// C_x_i[quadPt] += origvaluesD10[quadPt][j] * C_Array_i[globDOF];
				// C_y_i[quadPt] += origvaluesD01[quadPt][j] * C_Array_i[globDOF];
			}
		}

		// for ( int j = 0 ; j < N_BaseFunct; j++) rhs[j] = 0;
		for (int a = 0; a < N_S; a++)
		{ //"a" loop
		  // double* C_Array_a = U_Mode + a*lenMode;
			// Check
			double *Mode_Comp1_a = U_Mode + a * 2 * lenMode;	   // col Major
			double *Mode_Comp2_a = U_Mode + (a * 2 + 1) * lenMode; // col Major
			// double* phi_Array_a = Phi_Array + a*lenMode;??
			double *phi_Array_a = Phi_Array + a * lenPhi;

			double U1_Mode_a[N_Points2];
			double U1x_Mode_a[N_Points2];
			double U1y_Mode_a[N_Points2];
			double U1xx_Mode_a[N_Points2];
			double U1yy_Mode_a[N_Points2];

			double U2_Mode_a[N_Points2];
			double U2x_Mode_a[N_Points2];
			double U2y_Mode_a[N_Points2];
			double U2xx_Mode_a[N_Points2];
			double U2yy_Mode_a[N_Points2];

			for (int quadPt = 0; quadPt < N_Points2; quadPt++) // Initialize
			{
				// C_i[quadPt] = 0;
				U1_Mode_a[quadPt] = 0;
				U1x_Mode_a[quadPt] = 0;
				U1y_Mode_a[quadPt] = 0;
				U1xx_Mode_a[quadPt] = 0;
				U1yy_Mode_a[quadPt] = 0;

				U2_Mode_a[quadPt] = 0;
				U2x_Mode_a[quadPt] = 0;
				U2y_Mode_a[quadPt] = 0;
				U2xx_Mode_a[quadPt] = 0;
				U2yy_Mode_a[quadPt] = 0;
			}
			// for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++) C_x_a[quadPt] = 0;
			// for ( int quadPt = 0 ; quadPt < N_Points2; quadPt++) C_y_a[quadPt] = 0;

			// Obtain all values for C_a
			for (int quadPt = 0; quadPt < N_Points2; quadPt++)
			{
				for (int j = 0; j < N_BaseFunct; j++)
				{
					int globDOF = DOF[j];
					U1_Mode_a[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp1_a[globDOF];
					U1x_Mode_a[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp1_a[globDOF];
					U1y_Mode_a[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp1_a[globDOF];
					U1xx_Mode_a[quadPt] += origvaluesD20[quadPt][j] * Mode_Comp1_a[globDOF];
					U1yy_Mode_a[quadPt] += origvaluesD02[quadPt][j] * Mode_Comp1_a[globDOF];

					U2_Mode_a[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp2_a[globDOF];
					U2x_Mode_a[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp2_a[globDOF];
					U2y_Mode_a[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp2_a[globDOF];
					U2xx_Mode_a[quadPt] += origvaluesD20[quadPt][j] * Mode_Comp2_a[globDOF];
					U2yy_Mode_a[quadPt] += origvaluesD02[quadPt][j] * Mode_Comp2_a[globDOF];
				}
			}

			// Get Coefficients b1 and b2
			double *Param[MaxN_QuadPoints_2D];
			double **Coeffs = new double *[MaxN_QuadPoints_2D];
			for (int i = 0; i < MaxN_QuadPoints_2D; i++)
			{
				Coeffs[i] = new double[10]();
			}

			DO_Mode_Equation_Coefficients(N_Points2, X, Y, Param, Coeffs);

			for (int b = 0; b < N_S; b++)
			{ //"b" loop
				val = 0;
				double *Mode_Comp1_b = U_Mode + b * 2 * lenMode;	   // col Major
				double *Mode_Comp2_b = U_Mode + (b * 2 + 1) * lenMode; // col Major
				double *phi_Array_b = Phi_Array + b * lenPhi;

				double U1_Mode_b[N_Points2];
				double U1x_Mode_b[N_Points2];
				double U1y_Mode_b[N_Points2];

				double U2_Mode_b[N_Points2];
				double U2x_Mode_b[N_Points2];
				double U2y_Mode_b[N_Points2];

				for (int quadPt = 0; quadPt < N_Points2; quadPt++) // Initialize
				{
					// C_i[quadPt] = 0;
					U1_Mode_b[quadPt] = 0;
					U1x_Mode_b[quadPt] = 0;
					U1y_Mode_b[quadPt] = 0;

					U2_Mode_b[quadPt] = 0;
					U2x_Mode_b[quadPt] = 0;
					U2y_Mode_b[quadPt] = 0;
				}

				for (int quadPt = 0; quadPt < N_Points2; quadPt++)
				{
					for (int j = 0; j < N_BaseFunct; j++)
					{
						int globDOF = DOF[j];
						U1_Mode_b[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp1_b[globDOF];
						U1x_Mode_b[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp1_b[globDOF];
						U1y_Mode_b[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp1_b[globDOF];

						U2_Mode_b[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp2_b[globDOF];
						U2x_Mode_b[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp2_b[globDOF];
						U2y_Mode_b[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp2_b[globDOF];
					}
				}

				for (int qdpt = 0; qdpt < N_Points2; qdpt++)
				{
					double Mult = Weights2[qdpt] * AbsDetjk[qdpt];

					val += (U1_Mode_a[qdpt] * U1x_Mode_b[qdpt] + U2_Mode_a[qdpt] * U1y_Mode_b[qdpt]) * U1_Mode_i[qdpt] * Mult;
					val += (U1_Mode_a[qdpt] * U2x_Mode_b[qdpt] + U2_Mode_a[qdpt] * U2y_Mode_b[qdpt]) * U2_Mode_i[qdpt] * Mult;
				}

				for (int r = 0; r < lenPhi; r++)
				{
					phi_New[r] += -1.0 * (phi_Array_a[r] * phi_Array_b[r] - TDatabase::ParamDB->COVARIANCE_MATRIX_DO[a * N_S + b]) * val;
				}
			} //"b" loop ends
			val = 0;
			for (int qdpt = 0; qdpt < N_Points2; qdpt++)
			{
				double nu = Coeffs[qdpt][0];
				double Mult = Weights2[qdpt] * AbsDetjk[qdpt];
				val += (U1_Mode_a[qdpt] * U1x_Mean[qdpt] + U2_Mode_a[qdpt] * U1y_Mean[qdpt] + U1_Mean[qdpt] * U1x_Mode_a[qdpt] + U2_Mean[qdpt] * U1y_Mode_a[qdpt] + nu * U1xx_Mode_a[qdpt] + nu * U1yy_Mode_a[qdpt]) * U1_Mode_i[qdpt] * Mult;

				val += (U1_Mode_a[qdpt] * U2x_Mean[qdpt] + U2_Mode_a[qdpt] * U2y_Mean[qdpt] + U1_Mean[qdpt] * U2x_Mode_a[qdpt] + U2_Mean[qdpt] * U2y_Mode_a[qdpt] + nu * U2xx_Mode_a[qdpt] + nu * U2yy_Mode_a[qdpt]) * U2_Mode_i[qdpt] * Mult;
			}

			for (int i = 0; i < lenPhi; i++)
			{
				phi_New[i] += val * phi_Array_a[i] * -1.0;
			}

		} //"a" loop ends

	} // cell loop

	double timeStep = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
	for (int i = 0; i < lenPhi; i++)
	{
		phi_Array_i[i] = phi_Old_i[i] + timeStep * phi_New[i];
	}

	delete[] phi_New;
}

void DO_Mode_RHS_Aux_Param(double *in, double *out)
{
	out = in;
}
