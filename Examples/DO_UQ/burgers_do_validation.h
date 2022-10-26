/**
 * @file burgers_do_test.h
 * @brief Purpose:     Example file for solving the set of dynamically orthogonal field
                       equations for Burgers' equation.
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
//  Dynamically Orthogonal Field Equation Solution of Burgers' Equation //
// ===========================================================================//

// =======================================================================
//
// Purpose:     Example file for solving the set of dynamically orthogonal field
//              equations for Burgers' equation.
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

#include <MacroCell.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>
#include <memory>

extern "C"
{
#include <gridgen.h>

    void triangulate(char *, struct triangulateio *,
                     struct triangulateio *, struct triangulateio *);
}

// Include files related to the cell looping for the RHS filling part.

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

#include <QuadAffin.h>
#include <QuadBilinear.h>

//========================================================================//
// 								Example File							  //
// =======================================================================//

#define __SIN3__

void ExampleFile()
{

    OutPut("Example: Burgers' Equation - Dynamically Orthogonal Field Equation Solution" << endl);
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1Mean(double x, double y, double *values)
{
    double t = TDatabase::TimeDB->CURRENTTIME;
    values[0] = 0;

    // values[0] = exp(-t)*sin(Pi*x);
}

void InitialU2Mean(double x, double y, double *values)
{
    double t = TDatabase::TimeDB->CURRENTTIME;
    values[0] = 0;

    // values[0] = exp(-t)*-Pi*y*cos(Pi*x);
}

void InitialU1(double x, double y, double *values)
{
    double t = TDatabase::TimeDB->CURRENTTIME;

    // values[0] = exp(-t)*sin(Pi*x);
    values[0] = 0;
}

void InitialU2(double x, double y, double *values)
{
    double t = TDatabase::TimeDB->CURRENTTIME;
    values[0] = 0;

    // values[0] = exp(-t)*-Pi*y*cos(Pi*x);
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
    double t = TDatabase::TimeDB->CURRENTTIME, t1 = exp(-t);

    values[0] = sin(Pi * x) * t1;
    values[1] = Pi * cos(Pi * x) * t1;
    values[2] = 0;
    values[3] = -Pi * Pi * sin(Pi * x) * t1;
}

void ExactU2(double x, double y, double *values)
{
    double t = TDatabase::TimeDB->CURRENTTIME, t1 = exp(-t);

    values[0] = -Pi * y * cos(Pi * x) * t1;
    values[1] = Pi * Pi * y * sin(Pi * x) * t1;
    values[2] = -Pi * cos(Pi * x) * t1;
    values[3] = Pi * Pi * Pi * y * cos(Pi * x) * t1;
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
    double t = TDatabase::TimeDB->CURRENTTIME, t1 = exp(-t);
    value = 0;
    // switch(BdComp)
    // {
    //   case 0: value=sin(Pi*Param)*t1;
    //           break;
    //   case 1: value=0;
    //           break;
    //   case 2: value=sin(Pi*(1-Param))*t1;
    //           break;
    //   case 3: value=0;
    //           break;
    //   default: cout << "wrong boundary part number" << endl;
    //           break;
    // }
    // return;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
    double t = TDatabase::TimeDB->CURRENTTIME, t1 = exp(-t);
    value = 0;
    // switch(BdComp)
    // {
    //   case 0: value=0;
    //           break;
    //   case 1: value=Pi*Param*t1;
    //           break;
    //   case 2: value=-Pi*cos(Pi*(1-Param))*t1;
    //           break;
    //   case 3: value=-Pi*(1-Param)*t1;
    //           break;
    //   default: cout << "wrong boundary part number" << endl;
    //           break;
    // }
    // return;
}

// data on each quadrature point for mean
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
    double nu = 1. / TDatabase::ParamDB->RE_NR;
    int i;
    double *coeff, x, y;
    double u1, u1x, u1y, u2, u2x, u2y, u1t, u2t;
    double u1lap, u2lap, px, py;
    double t = TDatabase::TimeDB->CURRENTTIME, t1 = exp(-t), t2 = -exp(-t);

    for (i = 0; i < n_points; i++)
    {
        coeff = coeffs[i];

        x = X[i];
        y = Y[i];

        // prescribed solution
        u1 = sin(Pi * x) * t1;
        u1t = sin(Pi * x) * t2;
        u1x = Pi * cos(Pi * x) * t1;
        u1y = 0;
        u1lap = -Pi * Pi * sin(Pi * x) * t1;

        u2 = -Pi * y * cos(Pi * x) * t1;
        u2t = -Pi * y * cos(Pi * x) * t2;
        u2x = Pi * Pi * y * sin(Pi * x) * t1;
        u2y = -Pi * cos(Pi * x) * t1;
        u2lap = Pi * Pi * Pi * y * cos(Pi * x) * t1;

        coeff[0] = nu; // D(u):D(v) term
        // coeff[1] = -nu*u1lap  + u1t + u1*u1x + u2*u1y;
        // coeff[2] = -nu*u2lap  + u2t + u1*u2x + u2*u2y;
        coeff[1] = 0;
        coeff[2] = 0;
    }
}

void DO_Mean_Equation_Coefficients(int n_points, double *X, double *Y,
                                   double **parameters, double **coeffs)
{
    double nu = 1. / TDatabase::ParamDB->RE_NR;
    int i;
    double *coeff, x, y;
    double u1, u1x, u1y, u2, u2x, u2y, u1t, u2t;
    double u1lap, u2lap, px, py;
    double t = TDatabase::TimeDB->CURRENTTIME, t1 = exp(-t), t2 = -exp(-t);

    for (i = 0; i < n_points; i++)
    {
        coeff = coeffs[i];

        x = X[i];
        y = Y[i];

        // prescribed solution
        u1 = sin(Pi * x) * t1;
        u1t = sin(Pi * x) * t2;
        u1x = Pi * cos(Pi * x) * t1;
        u1y = 0;
        u1lap = -Pi * Pi * sin(Pi * x) * t1;

        u2 = -Pi * y * cos(Pi * x) * t1;
        u2t = -Pi * y * cos(Pi * x) * t2;
        u2x = Pi * Pi * y * sin(Pi * x) * t1;
        u2y = -Pi * cos(Pi * x) * t1;
        u2lap = Pi * Pi * Pi * y * cos(Pi * x) * t1;

        coeff[0] = nu; // D(u):D(v) term
        // coeff[1] = -nu*u1lap  + u1t + u1*u1x + u2*u1y;
        // coeff[2] = -nu*u2lap  + u2t + u1*u2x + u2*u2y;
        coeff[1] = 0;
        coeff[2] = 0;
    }
}

void DO_Mode_Equation_Coefficients(int n_points, double *X, double *Y,
                                   double **parameters, double **coeffs)
{
    double nu = 1. / TDatabase::ParamDB->RE_NR;
    int i;
    double *coeff, x, y;
    double u1, u1x, u1y, u2, u2x, u2y, u1t, u2t;
    double u1lap, u2lap, px, py;
    double t = TDatabase::TimeDB->CURRENTTIME, t1 = exp(-t), t2 = -exp(-t);

    for (i = 0; i < n_points; i++)
    {
        coeff = coeffs[i];

        x = X[i];
        y = Y[i];

        // prescribed solution
        u1 = sin(Pi * x) * t1;
        u1t = sin(Pi * x) * t2;
        u1x = Pi * cos(Pi * x) * t1;
        u1y = 0;
        u1lap = -Pi * Pi * sin(Pi * x) * t1;

        u2 = -Pi * y * cos(Pi * x) * t1;
        u2t = -Pi * y * cos(Pi * x) * t2;
        u2x = Pi * Pi * y * sin(Pi * x) * t1;
        u2y = -Pi * cos(Pi * x) * t1;
        u2lap = Pi * Pi * Pi * y * cos(Pi * x) * t1;

        coeff[0] = nu; // D(u):D(v) term
        // coeff[1] = -nu*u1lap  + u1t + u1*u1x + u2*u1y;
        // coeff[2] = -nu*u2lap  + u2t + u1*u2x + u2*u2y;
        coeff[1] = 0;
        coeff[2] = 0;
    }
}

// data on each quadrature point for mode
void LinCoeffs_Mode(int n_points, double *X, double *Y,
                    double **parameters, double **coeffs)
{
    double nu = 1. / TDatabase::ParamDB->RE_NR;
    int i;
    double *coeff, x, y;
    double u1, u1x, u1y, u2, u2x, u2y, u1t, u2t;
    double u1lap, u2lap, px, py;
    double t = TDatabase::TimeDB->CURRENTTIME, t1 = exp(-t), t2 = -exp(-t);

    for (i = 0; i < n_points; i++)
    {
        coeff = coeffs[i];

        x = X[i];
        y = Y[i];

        // prescribed solution
        u1 = sin(Pi * x) * t1;
        u1t = sin(Pi * x) * t2;
        u1x = Pi * cos(Pi * x) * t1;
        u1y = 0;
        u1lap = -Pi * Pi * sin(Pi * x) * t1;

        u2 = -Pi * y * cos(Pi * x) * t1;
        u2t = -Pi * y * cos(Pi * x) * t2;
        u2x = Pi * Pi * y * sin(Pi * x) * t1;
        u2y = -Pi * cos(Pi * x) * t1;
        u2lap = Pi * Pi * Pi * y * cos(Pi * x) * t1;

        coeff[0] = nu; // D(u):D(v) term
        coeff[1] = -nu * u1lap + u1t + u1 * u1x + u2 * u1y;
        coeff[2] = -nu * u2lap + u2t + u1 * u2x + u2 * u2y;
    }
}








// void normalizeStochasticModes(TFESpace2D *Fespace, TFEVectFunct2D *FeVector_Cmode, int N_S, double *C_Stoch_Norm)
// {

//     double *IPMatxFE = new double[N_S * N_S]();
//     int N_Cells = Fespace->GetN_Cells();
//     TCollection *coll = Fespace->GetCollection();

//     // Get the Global DOF arrays INdex from the FE Space.
//     int *GlobalNumbers = Fespace->GetGlobalNumbers();
//     int *BeginIndex = Fespace->GetBeginIndex();

//     // --- Quadrature Formula Arrays  ------------------//
//     int N_Points2;
//     double *Weights2, *t1, *t2; // Weights - Quadrature Weights array ; t1  -- Quadrature point ( xi ) in ref coordinate ; t2 -  Quadrature Point ( eta ) in Ref Coordinate
//     bool Needs2ndDer[1];
//     Needs2ndDer[0] = TRUE;
//     double AbsDetjk[MaxN_PointsForNodal2D];
//     double X[MaxN_PointsForNodal2D];
//     double Y[MaxN_PointsForNodal2D];

//     // FE Values Arrays
//     double **origvaluesD00; // Shape function values at quadrature Points
//     double **origvaluesD10; // Shape Function Derivatives ( x ) at Quadrature Points
//     double **origvaluesD01; // Shape Function Derivatives ( y ) at Quadrature Points
//     double **origvaluesD20; // Shape Function 2nd Derivatives ( x ) at Quadrature Points
//     double **origvaluesD02; // Shape Function 2nd Derivatives ( y ) at Quadrature Points

//     int lenMode = FeVector_Cmode->GetLength();
//     double *C_Mode_Array = FeVector_Cmode->GetValues();

//     double *C_Mode_Array_i = new double[lenMode]();
//     double *C_Mode_Array_j = new double[lenMode]();

//     double *stochFactors = new double[lenMode * N_S]();

//     for (int cellId = 0; cellId < N_Cells; cellId++)
//     { // cell loop
//         TBaseCell *currentCell = coll->GetCell(cellId);
//         // Get the "ID" of Finite Element for the given 2D Element ( Conforming/NonConforming-Order Finite Element : eg : it could be Conforming-2nd order Finite Element )
//         FE2D elementId = Fespace->GetFE2D(cellId, currentCell);
//         // Get the Class object for that 2d FEM Element , which has all the details like Shape functions , Nodal point locations for that location, Reference Transformation ( Affine , Bilinear )
//         TFE2D *element = TFEDatabase2D::GetFE2D(elementId);
//         TFEDesc2D *fedesc = element->GetFEDesc2D();
//         // Class for basis functions in 2D ( Number of basis functions ), basis function values and Derivatives
//         TBaseFunct2D *bf = element->GetBaseFunct2D();
//         // Get the Reference Elemet
//         BF2DRefElements RefElement = TFEDatabase2D::GetRefElementFromFE2D(elementId);
//         // Get the reference Transformation -- Affine Mapping / Bilnea Mapping of Triangle or Quadrilateral
//         RefTrans2D referenceTransformation = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(elementId);
//         // Get the number of basis functions in the Current Cell ( Number of Local DOF)
//         int N_BaseFunct = element->GetN_DOF();
//         // Type of Basis Function in 2D
//         BaseFunct2D BaseFunct_ID = element->GetBaseFunct2D_ID();

//         // get cell measure
//         double hK = currentCell->GetDiameter();

//         switch (referenceTransformation)
//         {
//         case QuadBilinear:
//         {
//             int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
//             QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
//             TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
//             QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

//             // Set the values on the Reference Cell
//             TRefTrans2D *F_K = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
//             TFEDatabase2D::SetCellForRefTrans(currentCell, QuadBilinear); // Set the Cell for Current reference Transformation

//             // Get Original Coordinates from reference Coordinates and the Determinant of jacobian
//             TFEDatabase2D::GetOrigFromRef(QuadBilinear, N_Points2, t1, t2, X, Y, AbsDetjk); // Get the Original Co-orinates for the cell from xi values

//             // Get all the original Values from the Referece cell values.
//             TFEDatabase2D::GetOrigValues(QuadBilinear, 1, &BaseFunct_ID, N_Points2, t1, t2, QF2, Needs2ndDer);

//             // The below are 2D arrays in the form
//             // Values[QuadraturePointLocation][ShapeFunction]  i.e, the Value of Shapefunction at all quadrature points for each shape functions
//             origvaluesD00 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D00); // Shape Function Values at Quadrature Points
//             origvaluesD10 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D10); // Shape Function Derivative Values at Quadrature Points
//             origvaluesD01 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D01); // Shape Function Derivative Values at Quadrature Point
//             origvaluesD20 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D20); // Shape Function 2nd Derivative Values at Quadrature Point
//             origvaluesD02 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D02); // Shape Function 2nd Derivative Values at Quadrature Point
//             break;
//         }

//         case QuadAffin:
//         {
//             int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
//             QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
//             TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
//             QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

//             // Set the values on the Reference Cell
//             TRefTrans2D *F_K = TFEDatabase2D::GetRefTrans2D(QuadAffin);
//             TFEDatabase2D::SetCellForRefTrans(currentCell, QuadAffin); // Set the Cell for Current reference Transformation

//             // Get Original Coordinates from reference Coordinates and the Determinant of jacobian
//             TFEDatabase2D::GetOrigFromRef(QuadAffin, N_Points2, t1, t2, X, Y, AbsDetjk); // Get the Original Co-orinates for the cell from xi values

//             // Get all the original Values from the Referece cell values.
//             TFEDatabase2D::GetOrigValues(QuadAffin, 1, &BaseFunct_ID, N_Points2, t1, t2, QF2, Needs2ndDer);

//             // The below are 2D arrays in the form
//             // Values[QuadraturePointLocation][ShapeFunction]  i.e, the Value of Shapefunction at all quadrature points for each shape functions
//             origvaluesD00 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D00); // Shape Function Values at Quadrature Points
//             origvaluesD10 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D10); // Shape Function Derivative Values at Quadrature Points
//             origvaluesD01 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D01); // Shape Function Derivative Values at Quadrature Point
//             origvaluesD20 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D20); // Shape Function 2nd Derivative Values at Quadrature Points
//             origvaluesD02 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D02); // Shape Function 2nd Derivative Values at Quadrature Point

//             break;
//         }

//         default:
//         {
//             cout << " [ERROR] - Error in File : CoeffEqn_DO.C " << endl;
//             cout << " Unknown Reftype " << endl;
//             cout << " REF TYPE : " << referenceTransformation << endl;
//             exit(0);
//             break;
//         }
//         }

//         int *DOF = GlobalNumbers + BeginIndex[cellId];

//         for (int i = 0; i < N_S; i++)
//         { // i-loop start
//             memcpy(C_Mode_Array_i, C_Mode_Array + (2 * i * lenMode), lenMode * SizeOfDouble);
//             double C_i[N_Points2];

//             for (int quadPt = 0; quadPt < N_Points2; quadPt++)
//                 C_i[quadPt] = 0.;
//             // Obtain all values for C_a

//             for (int quadPt = 0; quadPt < N_Points2; quadPt++)
//             { // quadpt loop start
//                 for (int nb = 0; nb < N_BaseFunct; nb++)
//                 {
//                     int globDOF = DOF[nb];
//                     C_i[quadPt] += origvaluesD00[quadPt][nb] * C_Mode_Array_i[globDOF];
//                     stochFactors[globDOF + i * lenMode] += 1;
//                 }
//             } // quadpt loop end

//             for (int j = 0; j < N_S; j++)
//             { // j-loop start

//                 memcpy(C_Mode_Array_j, C_Mode_Array + (2 * j * lenMode), lenMode * SizeOfDouble);
//                 double C_j[N_Points2];
//                 for (int quadPt = 0; quadPt < N_Points2; quadPt++)
//                     C_j[quadPt] = 0.;
//                 // Obtain all values for C_a
//                 for (int quadPt = 0; quadPt < N_Points2; quadPt++)
//                 { // quadpt loop start
//                     for (int nb = 0; nb < N_BaseFunct; nb++)
//                     {
//                         int globDOF = DOF[nb];
//                         C_j[quadPt] += origvaluesD00[quadPt][nb] * C_Mode_Array_j[globDOF];
//                     }
//                 } // quadpt loop end

//                 // INner Quadrature Loop
//                 for (int quadPt = 0; quadPt < N_Points2; quadPt++)
//                 {
//                     double Mult = Weights2[quadPt] * AbsDetjk[quadPt];
//                     IPMatxFE[i + j * N_S] += C_i[quadPt] * C_j[quadPt] * Mult;

//                 } // Inner Quadrature Loop

//             } // j-loop end

//         } // i-loop end

//         // --
//     } // cell loop

//     for (int j = 0; j < N_S; j++)
//     {
//         for (int i = 0; i < lenMode; i++)
//         {
//             C_Stoch_Norm[i + (2 * j * lenMode)] = C_Mode_Array[i + (2 * j * lenMode)] / sqrt(IPMatxFE[j + j * N_S]);
//         }
//     }

//     double *C_temp = new double[lenMode*N_S]();
//      for (int j = 0; j < N_S; j++)
//     {
//         for (int i = 0; i < lenMode; i++)
//         {
//             C_temp[i + j * lenMode] = C_Stoch_Norm[i + (2 * j * lenMode)];
//         }
//     }
//     double *IPNew = new double[N_S * N_S]();
//     calcIPMatx(IPNew, C_temp, lenMode, N_S, 'C');
//     for (int j = 0; j < N_S; j++)
//     {
//         for (int i = 0; i < lenMode; i++)
//         {
//             C_Stoch_Norm[i + (2 * j * lenMode)] = C_temp[i+j*lenMode]/ sqrt(IPNew[j + j * N_S]);
//         }
//     }

//     delete[] C_Mode_Array_i;
//     delete[] C_Mode_Array_j;

//     delete[] IPMatxFE;

//     return;
// } // stochastic normalization function end


void normalizeStochasticMean(TFESpace2D *Fespace, TFEVectFunct2D *FeVector_CMean, int N_S, double *C_Stoch_Norm)
{
    double *IPMatxFE = new double[N_S * N_S]();
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

    int lenMean = FeVector_CMean->GetLength();
    double *C_Mean_Array = FeVector_CMean->GetValues();

    double *C_Mean = new double[lenMean]();

    double mfeval = 0.0;

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
            int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
            QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
            TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
            QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

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
            int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
            QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
            TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
            QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

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

        memcpy(C_Mean, C_Mean_Array, lenMean * SizeOfDouble);
        double C_Mean_1[N_Points2];

        for (int quadPt = 0; quadPt < N_Points2; quadPt++)
            C_Mean_1[quadPt] = 0.;
        // Obtain all values for C_a

        for (int quadPt = 0; quadPt < N_Points2; quadPt++)
        { // quadpt loop start
            for (int nb = 0; nb < N_BaseFunct; nb++)
            {
                int globDOF = DOF[nb];
                C_Mean_1[quadPt] += origvaluesD00[quadPt][nb] * C_Mean[globDOF];
            }
        } // quadpt loop end

        // INner Quadrature Loop
        for (int quadPt = 0; quadPt < N_Points2; quadPt++)
        {
            double Mult = Weights2[quadPt] * AbsDetjk[quadPt];
            mfeval += C_Mean_1[quadPt] * C_Mean_1[quadPt] * Mult;

        } // Inner Quadrature Loop

        // --
    } // cell loop

    double IPVal = 0.0;

    for (int i = 0; i < lenMean; i++)
    {
        IPVal += C_Mean[i] * C_Mean[i];
    }
    for (int i = 0; i < lenMean; i++)
    {
        C_Stoch_Norm[i] = C_Mean[i] * (sqrt(IPVal / mfeval));
    }
}

void DO_Mean_RHS(TFESpace2D *Fespace, TFEVectFunct2D *FeVector_Mode, int N_S, double *GlobalRhs_mean, int N_U)
{

    int N_Cells = Fespace->GetN_Cells();
    TCollection *coll = Fespace->GetCollection();

    double *U_Mode = FeVector_Mode->GetValues();
    int lenMode = FeVector_Mode->GetLength();
    // cout << "Length Mode in DO Mean RHS: " << lenMode << endl;
    int lenMean = N_U;
    double *Mode_Comp1_a = new double[lenMode]();
    double *Mode_Comp2_a = new double[lenMode]();
    double *Mode_Comp1_b = new double[lenMode]();
    double *Mode_Comp2_b = new double[lenMode]();
    double val1 = 0;
    double val2 = 0;

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
    double ip = 0.0;
    for (int cellId = 0; cellId < N_Cells; cellId++)
    { // Cell Loop
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
            int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
            QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
            TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
            QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

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
            int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
            QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
            TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
            QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

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
        {                                                                                   //"a" loop
            memcpy(Mode_Comp1_a, U_Mode + (a * 2 * lenMode), lenMode * SizeOfDouble);       // col Major
            memcpy(Mode_Comp2_a, U_Mode + ((a * 2 + 1) * lenMode), lenMode * SizeOfDouble); // col Major
                                                                                            //  double* phi_Array_a = Phi_Array + a*lenMode;??

            double U1_Mode_a[N_Points2];
            double U2_Mode_a[N_Points2];

            for (int quadPt = 0; quadPt < N_Points2; quadPt++) // Initialize
            {
                // C_i[quadPt] = 0;
                U1_Mode_a[quadPt] = 0;
                U2_Mode_a[quadPt] = 0;
            }

            // Obtain all values for C_a
            for (int quadPt = 0; quadPt < N_Points2; quadPt++)
            {
                for (int j = 0; j < N_BaseFunct; j++)
                {
                    int globDOF = DOF[j];
                    U1_Mode_a[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp1_a[globDOF];
                    U2_Mode_a[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp2_a[globDOF];
                }
            }

            for (int b = 0; b < N_S; b++)
            { //"b" loop
                val1 = 0;
                val2 = 0;
                memcpy(Mode_Comp1_b, U_Mode + (b * 2 * lenMode), lenMode * SizeOfDouble);       // col Major
                memcpy(Mode_Comp2_b, U_Mode + ((b * 2 + 1) * lenMode), lenMode * SizeOfDouble); // col Major

                double U1x_Mode_b[N_Points2];
                double U1y_Mode_b[N_Points2];
                double U2x_Mode_b[N_Points2];
                double U2y_Mode_b[N_Points2];

                for (int quadPt = 0; quadPt < N_Points2; quadPt++) // Initialize
                {
                    // C_i[quadPt] = 0;
                    U1x_Mode_b[quadPt] = 0;
                    U1y_Mode_b[quadPt] = 0;
                    U2x_Mode_b[quadPt] = 0;
                    U2y_Mode_b[quadPt] = 0;
                }

                for (int quadPt = 0; quadPt < N_Points2; quadPt++)
                {
                    for (int j = 0; j < N_BaseFunct; j++)
                    {
                        int globDOF = DOF[j];
                        U1x_Mode_b[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp1_b[globDOF];
                        U1y_Mode_b[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp1_b[globDOF];
                        U2x_Mode_b[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp2_b[globDOF];
                        U2y_Mode_b[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp2_b[globDOF];
                    }
                }

                for (int qdpt = 0; qdpt < N_Points2; qdpt++)
                {
                    double Mult = Weights2[qdpt] * AbsDetjk[qdpt];
                    double *orgD00 = origvaluesD00[qdpt];

                    val1 = -1.0 * ((U1_Mode_a[qdpt] * U1x_Mode_b[qdpt]) + (U2_Mode_a[qdpt] * U1y_Mode_b[qdpt])) * TDatabase::ParamDB->COVARIANCE_MATRIX_DO[a * N_S + b] * Mult;

                    val2 = -1.0 * ((U1_Mode_a[qdpt] * U2x_Mode_b[qdpt]) + (U2_Mode_a[qdpt] * U2y_Mode_b[qdpt])) * TDatabase::ParamDB->COVARIANCE_MATRIX_DO[a * N_S + b] * Mult; // check if = or +=
                    
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
            GlobalRhs_mean[GlobalDOF + lenMean] += rhs2[j];
        }
        // --
        for (int i = 0; i < MaxN_QuadPoints_2D; i++)
        {
            delete[] Coeffs[i];
        }
        delete[] Coeffs;
    } // cell loop
    cout << "value of ip" << ip << endl;

    delete[] Mode_Comp1_a;
    delete[] Mode_Comp2_a;
    delete[] Mode_Comp1_b;
    delete[] Mode_Comp2_b;
} // looks right

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

    double *Mean_Comp1 = new double[lenMean]();
    memcpy(Mean_Comp1, U_Mean, lenMean * SizeOfDouble);
    double *Mean_Comp2 = new double[lenMean]();
    memcpy(Mean_Comp2, U_Mean + lenMean, lenMean * SizeOfDouble); //

    double *Mode_Comp1_i = new double[lenMode]();
    memcpy(Mode_Comp1_i, U_Mode + i_index * 2 * lenMode, lenMode * SizeOfDouble); // col Major
    double *Mode_Comp2_i = new double[lenMode]();
    memcpy(Mode_Comp2_i, U_Mode + (i_index * 2 + 1) * lenMode, lenMode * SizeOfDouble); // col Major

    double val1 = 0;
    double val2 = 0;
    double *Mode_Comp1_a = new double[lenMode]();
    double *Mode_Comp2_a = new double[lenMode]();
    double *Mode_Comp1_b = new double[lenMode]();
    double *Mode_Comp2_b = new double[lenMode]();
    double *Mode_Comp1_p = new double[lenMode]();
    double *Mode_Comp2_p = new double[lenMode]();

    double *ipval02 = new double[N_S]();
    double *ipval01 = new double[N_S]();

    double *ipval11 = new double[N_S]();
    double *ipval12 = new double[N_S]();
    for (int p = 0; p < N_S; p++)
    {
        memcpy(Mode_Comp1_p, U_Mode + (p * 2 * lenMode), lenMode * SizeOfDouble);
        memcpy(Mode_Comp2_p, U_Mode + ((p * 2 + 1) * lenMode), lenMode * SizeOfDouble);
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
                int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
                QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
                TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
                QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

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
                int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
                QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
                TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
                QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

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
            double *Param[MaxN_QuadPoints_2D];
            double **Coeffs = new double *[MaxN_QuadPoints_2D];
            for (int i = 0; i < MaxN_QuadPoints_2D; i++)
            {
                Coeffs[i] = new double[10]();
            }
            DO_Mode_Equation_Coefficients(N_Points2, X, Y, Param, Coeffs);
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

            double U1_Mode_p[N_Points2];
            double U2_Mode_p[N_Points2];

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

                U1_Mode_p[quadPt] = 0;
                U2_Mode_p[quadPt] = 0;
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

                    U1_Mode_p[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp1_p[globDOF];
                    U2_Mode_p[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp2_p[globDOF];
                }
            }
            for (int qdpt = 0; qdpt < N_Points2; qdpt++)
            {
                double Mult = Weights2[qdpt] * AbsDetjk[qdpt];
                double nu = Coeffs[qdpt][0]; // nu

                ipval01[p] += (-1.0 * (U1_Mode_i[qdpt] * U1x_Mean[qdpt] + U2_Mode_i[qdpt] * U1y_Mean[qdpt])) * U1_Mode_p[qdpt] * Mult;

                ipval01[p] += (-1.0 * (U1_Mean[qdpt] * U1x_Mode_i[qdpt] + U2_Mean[qdpt] * U1y_Mode_i[qdpt])) * U1_Mode_p[qdpt] * Mult;

                ipval01[p] += nu * (U1xx_Mode_i[qdpt] + U1yy_Mode_i[qdpt]) * U1_Mode_p[qdpt] * Mult;

                ipval02[p] += (-1.0 * (U1_Mode_i[qdpt] * U2x_Mean[qdpt] + U2_Mode_i[qdpt] * U2y_Mean[qdpt])) * U2_Mode_p[qdpt] * Mult;

                ipval02[p] += (-1.0 * (U1_Mean[qdpt] * U2x_Mode_i[qdpt] + U2_Mean[qdpt] * U2y_Mode_i[qdpt])) * U2_Mode_p[qdpt] * Mult;

                ipval02[p] += nu * (U2xx_Mode_i[qdpt] + U2yy_Mode_i[qdpt]) * U2_Mode_p[qdpt] * Mult;
            }
            for (int i = 0; i < MaxN_QuadPoints_2D; i++)
            {
                delete[] Coeffs[i];
            }
            delete[] Coeffs;
        } // cell loop end


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
                int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
                QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
                TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
                QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

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
                int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
                QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
                TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
                QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

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
            double *Param[MaxN_QuadPoints_2D];
            double **Coeffs = new double *[MaxN_QuadPoints_2D];
            for (int i = 0; i < MaxN_QuadPoints_2D; i++)
            {
                Coeffs[i] = new double[10]();
            }
            DO_Mode_Equation_Coefficients(N_Points2, X, Y, Param, Coeffs);
            double U1_Mode_p[N_Points2];
            double U2_Mode_p[N_Points2];
            for (int quadPt = 0; quadPt < N_Points2; quadPt++) // Initialize
            {
                U1_Mode_p[N_Points2] = 0;
                U2_Mode_p[N_Points2] = 0;
            }
            for (int quadPt = 0; quadPt < N_Points2; quadPt++)
            {
                for (int j = 0; j < N_BaseFunct; j++)
                {
                    int globDOF = DOF[j];
                    U1_Mode_p[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp1_p[globDOF];
                    U2_Mode_p[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp2_p[globDOF];
                }
            }
            for (int a = 0; a < N_S; a++)
            {                                                                                 // a loop for ip calc 2
                memcpy(Mode_Comp1_a, U_Mode + a * 2 * lenMode, lenMode * SizeOfDouble);       // col Major
                memcpy(Mode_Comp2_a, U_Mode + (a * 2 + 1) * lenMode, lenMode * SizeOfDouble); // col Major

                double U1_Mode_a[N_Points2];
                double U2_Mode_a[N_Points2];

                for (int quadPt = 0; quadPt < N_Points2; quadPt++) // Initialize
                {
                    U1_Mode_a[quadPt] = 0;
                    U2_Mode_a[quadPt] = 0;
                }

                // Obtain all values for U_a
                for (int quadPt = 0; quadPt < N_Points2; quadPt++)
                {
                    for (int j = 0; j < N_BaseFunct; j++)
                    {
                        int globDOF = DOF[j];
                        U1_Mode_a[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp1_a[globDOF];
                        U2_Mode_a[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp2_a[globDOF];
                    }
                }

                for (int b = 0; b < N_S; b++)
                {                                                                                 // b loop start for ip val calc 2
                    memcpy(Mode_Comp1_b, U_Mode + b * 2 * lenMode, lenMode * SizeOfDouble);       // col Major
                    memcpy(Mode_Comp2_b, U_Mode + (b * 2 + 1) * lenMode, lenMode * SizeOfDouble); // col Major

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

                    for (int c = 0; c < N_S; c++)
                    {
                        for (int qdpt = 0; qdpt < N_Points2; qdpt++)
                        {
                            double Mult = Weights2[qdpt] * AbsDetjk[qdpt];

                            ipval11[p] += -1.0 * (TDatabase::ParamDB->COVARIANCE_INVERSE_DO[N_S * c + i_index] * TDatabase::ParamDB->COSKEWNESS_MATRIX_DO[N_S * N_S * b + N_S * c + a]) * (U1_Mode_a[qdpt] * U1x_Mode_b[qdpt] + U2_Mode_a[qdpt] * U1y_Mode_b[qdpt]) * U1_Mode_p[qdpt] * Mult;

                            ipval12[p] += -1.0 * (TDatabase::ParamDB->COVARIANCE_INVERSE_DO[N_S * c + i_index] * TDatabase::ParamDB->COSKEWNESS_MATRIX_DO[N_S * N_S * b + N_S * c + a]) * (U1_Mode_a[qdpt] * U2x_Mode_b[qdpt] + U2_Mode_a[qdpt] * U2y_Mode_b[qdpt]) * U2_Mode_p[qdpt] * Mult;
                        }
                    }//c loop end for ipval cal2

                } // b loop end for ipval calc 2
            }     // a loop end for ip calc 2
            for (int i = 0; i < MaxN_QuadPoints_2D; i++)
            {
                delete[] Coeffs[i];
            }
            delete[] Coeffs;
        } // cell loop

    } // p loop end
    
    for (int cellId = 0; cellId < N_Cells; cellId++)
    { // cell loop final rhs ass 1
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
            int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
            QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
            TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
            QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

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
            int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
            QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
            TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
            QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

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
        DO_Mode_Equation_Coefficients(N_Points2, X, Y, Param, Coeffs);

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

        double rhs1[N_BaseFunct];
        double rhs2[N_BaseFunct];
        for (int j = 0; j < N_BaseFunct; j++)
        {
            rhs1[j] = 0;
            rhs2[j] = 0;
        }

        for (int qdpt = 0; qdpt < N_Points2; qdpt++)
        {
            double Mult = Weights2[qdpt] * AbsDetjk[qdpt];
            double *orgD00 = origvaluesD00[qdpt];
            double *orgD10 = origvaluesD10[qdpt];
            double *orgD01 = origvaluesD01[qdpt];
            double *orgD20 = origvaluesD20[qdpt];
            double *orgD02 = origvaluesD02[qdpt];
            double nu = Coeffs[qdpt][0]; // nu

            val1 = -1.0 * (U1_Mode_i[qdpt] * U1x_Mean[qdpt] + U2_Mode_i[qdpt] * U1y_Mean[qdpt]) * Mult;
            val1 += -1.0 * (U1_Mean[qdpt] * U1x_Mode_i[qdpt] + U2_Mode_i[qdpt] * U1y_Mode_i[qdpt]) * Mult;
            val1 += (U1_Mode_i[qdpt] * U1x_Mode_i[qdpt] + U2_Mode_i[qdpt] * U1y_Mode_i[qdpt]) * Mult;
            val1 += nu * (U1xx_Mode_i[qdpt] + U1yy_Mode_i[qdpt]) * Mult;

            val2 = -1.0 * (U1_Mode_i[qdpt] * U2x_Mean[qdpt] + U2_Mode_i[qdpt] * U2y_Mean[qdpt]) * Mult;
            val2 += -1.0 * (U1_Mean[qdpt] * U2x_Mode_i[qdpt] + U2_Mode_i[qdpt] * U2y_Mode_i[qdpt]) * Mult;
            val2 += (U1_Mode_i[qdpt] * U2x_Mode_i[qdpt] + U2_Mode_i[qdpt] * U2y_Mode_i[qdpt]) * Mult;
            val2 += nu * (U2xx_Mode_i[qdpt] + U2yy_Mode_i[qdpt]) * Mult;
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
        for (int i = 0; i < MaxN_QuadPoints_2D; i++)
        {
            delete[] Coeffs[i];
        }
        delete[] Coeffs;
    } // cell loop final rhs ass 1 end

    for (int cellId = 0; cellId < N_Cells; cellId++)
    { // cell loop final rhs ass 2
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
            int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
            QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
            TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
            QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

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
            int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
            QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
            TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
            QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

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
        DO_Mode_Equation_Coefficients(N_Points2, X, Y, Param, Coeffs);

        // Save Values of C at all quadrature points for I component

        double rhs1[N_BaseFunct];
        double rhs2[N_BaseFunct];
        for (int j = 0; j < N_BaseFunct; j++)
        {
            rhs1[j] = 0;
            rhs2[j] = 0;
        }

        for (int a = 0; a < N_S; a++)
        {
            memcpy(Mode_Comp1_a, U_Mode + a * 2 * lenMode, lenMode * SizeOfDouble);       // col Major
            memcpy(Mode_Comp2_a, U_Mode + (a * 2 + 1) * lenMode, lenMode * SizeOfDouble); // col Major

            double U1_Mode_a[N_Points2];
            double U2_Mode_a[N_Points2];

            for (int quadPt = 0; quadPt < N_Points2; quadPt++) // Initialize
            {
                // C_i[quadPt] = 0;
                U1_Mode_a[quadPt] = 0;
                U2_Mode_a[quadPt] = 0;
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
                    U2_Mode_a[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp2_a[globDOF];
                }
            }
            for (int b = 0; b < N_S; b++)
            {
                memcpy(Mode_Comp1_b, U_Mode + b * 2 * lenMode, SizeOfDouble);                 // col Major
                memcpy(Mode_Comp2_b, U_Mode + (b * 2 + 1) * lenMode, lenMode * SizeOfDouble); // col Major
                double U1x_Mode_b[N_Points2];
                double U1y_Mode_b[N_Points2];

                double U2x_Mode_b[N_Points2];
                double U2y_Mode_b[N_Points2];

                for (int quadPt = 0; quadPt < N_Points2; quadPt++) // Initialize
                {
                    U1x_Mode_b[quadPt] = 0;
                    U1y_Mode_b[quadPt] = 0;

                    U2x_Mode_b[quadPt] = 0;
                    U2y_Mode_b[quadPt] = 0;
                }

                for (int quadPt = 0; quadPt < N_Points2; quadPt++)
                {
                    for (int j = 0; j < N_BaseFunct; j++)
                    {
                        int globDOF = DOF[j];
                        U1x_Mode_b[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp1_b[globDOF];
                        U1y_Mode_b[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp1_b[globDOF];

                        U2x_Mode_b[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp2_b[globDOF];
                        U2y_Mode_b[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp2_b[globDOF];
                    }
                }

                for (int c = 0; c < N_S; c++)
                {
                    for (int qdpt = 0; qdpt < N_Points2; qdpt++)
                    {
                        double Mult = Weights2[qdpt] * AbsDetjk[qdpt];
                        double *orgD00 = origvaluesD00[qdpt];

                        val1 = -1.0 * (TDatabase::ParamDB->COVARIANCE_INVERSE_DO[N_S * c + i_index] * TDatabase::ParamDB->COSKEWNESS_MATRIX_DO[N_S * N_S * b + N_S * c + a]) * (U1_Mode_a[qdpt] * U1x_Mode_b[qdpt] + U2_Mode_a[qdpt] * U1y_Mode_b[qdpt]) * Mult;

                        val2 = -1.0 * (TDatabase::ParamDB->COVARIANCE_INVERSE_DO[N_S * c + i_index] * TDatabase::ParamDB->COSKEWNESS_MATRIX_DO[N_S * N_S * b + N_S * c + a]) * (U1_Mode_a[qdpt] * U2x_Mode_b[qdpt] + U2_Mode_a[qdpt] * U2y_Mode_b[qdpt]) * Mult;
                        for (int j = 0; j < N_BaseFunct; j++)
                        {
                            rhs1[j] += val1 * orgD00[j]; // * Mult;
                            rhs2[j] += val2 * orgD00[j]; // * Mult;
                        }
                    }

                    for (int j = 0; j < N_BaseFunct; j++)
                    {
                        int GlobalDOF = DOF[j];
                        GlobalRhs_mode[GlobalDOF] += rhs1[j];
                        GlobalRhs_mode[GlobalDOF + lenMode] += rhs2[j];
                    }
                } // cc loop end
            }     // b
        }         // a
        for (int i = 0; i < MaxN_QuadPoints_2D; i++)
        {
            delete[] Coeffs[i];
        }
        delete[] Coeffs;
    } // cell loop final rhs ass 2 end

    for (int cellId = 0; cellId < N_Cells; cellId++)
    { // cell loop rhs ass 3
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
            int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
            QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
            TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
            QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

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
            int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
            QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
            TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
            QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

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
        DO_Mode_Equation_Coefficients(N_Points2, X, Y, Param, Coeffs);
        for (int p = 0; p < N_S; p++)
        {
            memcpy(Mode_Comp1_p, U_Mode + p * 2 * lenMode, lenMode * SizeOfDouble);       // col Major
            memcpy(Mode_Comp2_p, U_Mode + (p * 2 + 1) * lenMode, lenMode * SizeOfDouble); // col Major

            double U1_Mode_p[N_Points2];
            double U2_Mode_p[N_Points2];

            for (int quadPt = 0; quadPt < N_Points2; quadPt++) // Initialize
            {
                U1_Mode_p[quadPt] = 0;
                U2_Mode_p[quadPt] = 0;
            }

            double rhs1[N_BaseFunct];
            double rhs2[N_BaseFunct];
            for (int j = 0; j < N_BaseFunct; j++)
            {
                rhs1[j] = 0;
                rhs2[j] = 0;
            }

            for (int quadPt = 0; quadPt < N_Points2; quadPt++)
            {
                double Mult = Weights2[quadPt] * AbsDetjk[quadPt];
                double *orgD00 = origvaluesD00[quadPt];

                double nu = Coeffs[quadPt][0]; // nu

                val1 = -1.0 * (ipval01[p] + ipval11[p]) * U1_Mode_p[quadPt] * Mult;
                val2 = -1.0 * (ipval02[p] + ipval12[p]) * U2_Mode_p[quadPt] * Mult;

                for (int j = 0; j < N_BaseFunct; j++)
                {
                    rhs1[j] += val1 * orgD00[j]; // * Mult;
                    rhs2[j] += val2 * orgD00[j]; // * Mult;
                }

            } // Inner Quadrature Loop for p
            for (int j = 0; j < N_BaseFunct; j++)
            {
                int GlobalDOF = DOF[j];
                GlobalRhs_mode[GlobalDOF] += rhs1[j];
                GlobalRhs_mode[GlobalDOF + lenMode] += rhs2[j];
            }
        } // p loop
        for (int i = 0; i < MaxN_QuadPoints_2D; i++)
        {
            delete[] Coeffs[i];
        }
        delete[] Coeffs;
    } // cell loop rhs ass 3 end

    delete[] Mean_Comp1;
    delete[] Mean_Comp2;
    delete[] Mode_Comp1_i;
    delete[] Mode_Comp2_i;
    delete[] Mode_Comp1_a;
    delete[] Mode_Comp2_a;
    delete[] Mode_Comp1_b;
    delete[] Mode_Comp2_b;
    delete[] Mode_Comp1_p;
    delete[] Mode_Comp2_p;
    delete[] ipval01;
    delete[] ipval02;
    delete[] ipval11;
    delete[] ipval12;
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

    double val = 0;
    double *U_Mode = FeVector_Mode->GetValues();
    int lenMode = FeVector_Mode->GetLength();

    // Check
    double *U_Mean = FeVector_Mean->GetValues();
    int lenMean = FeVector_Mean->GetLength();

    double *Mean_Comp1 = new double[lenMean]();
    memcpy(Mean_Comp1, U_Mean, lenMean * SizeOfDouble);
    double *Mean_Comp2 = new double[lenMean]();
    memcpy(Mean_Comp2, U_Mean + lenMean, lenMean * SizeOfDouble); //

    double *Phi_Array = FEVector_Phi->GetValues();
    int lenPhi = FEVector_Phi->GetLength();

    double *Phi_Old = new double[N_R * N_S]();
    memcpy(Phi_Old, Phi_Array, N_R * N_S * SizeOfDouble);

    double *phi_New = new double[lenPhi]();

    // Check
    double *Mode_Comp1_i = new double[lenMode]();
    memcpy(Mode_Comp1_i, U_Mode + i_index * 2 * lenMode, lenMode * SizeOfDouble); // col Major
    double *Mode_Comp2_i = new double[lenMode]();
    memcpy(Mode_Comp2_i, U_Mode + (i_index * 2 + 1) * lenMode, lenMode * SizeOfDouble); // col Major
    // 2nd Component??

    double *phi_Array_i = Phi_Array + i_index * lenPhi;
    double *phi_Old_i = new double[lenPhi]();
    memcpy(phi_Old_i, Phi_Old + (i_index * lenPhi), lenPhi * SizeOfDouble);

    //   double* phi_New = new double[lenPhi]();

    double *Mode_Comp1_a = new double[lenMode]();
    double *Mode_Comp2_a = new double[lenMode]();
    double *phi_Array_a = new double[lenPhi]();

    double *Mode_Comp1_b = new double[lenMode]();
    double *Mode_Comp2_b = new double[lenMode]();
    double *phi_Array_b = new double[lenPhi]();

    double ipval1 = 0.0;
    double ipval2 = 0.0;
    for (int a = 0; a < N_S; a++)
    { // a loop ip1
        for (int cellId = 0; cellId < N_Cells; cellId++)
        { // cell loop ip1
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
                int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
                QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
                TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
                QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

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
                int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
                QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
                TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
                QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

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

            // Get Coefficients b1 and b2
            double *Param[MaxN_QuadPoints_2D];
            double **Coeffs = new double *[MaxN_QuadPoints_2D];
            for (int i = 0; i < MaxN_QuadPoints_2D; i++)
            {
                Coeffs[i] = new double[10]();
            }
            DO_Mode_Equation_Coefficients(N_Points2, X, Y, Param, Coeffs);

            memcpy(Mode_Comp1_a, U_Mode + (a * 2 * lenMode), lenMode * SizeOfDouble);
            memcpy(Mode_Comp2_a, U_Mode + ((a * 2 + 1) * lenMode), lenMode * SizeOfDouble);
            memcpy(phi_Array_a, Phi_Array + (a * lenPhi), lenPhi * SizeOfDouble);
            // Save Values of C at all quadrature points for I component
            double U1_Mode_i[N_Points2];
            double U2_Mode_i[N_Points2];

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

                U2_Mode_i[quadPt] = 0;

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

                    U2_Mode_i[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp2_i[globDOF];

                    U1_Mean[quadPt] += origvaluesD00[quadPt][j] * Mean_Comp1[globDOF];
                    U1x_Mean[quadPt] += origvaluesD10[quadPt][j] * Mean_Comp1[globDOF];
                    U1y_Mean[quadPt] += origvaluesD01[quadPt][j] * Mean_Comp1[globDOF];

                    U2_Mean[quadPt] += origvaluesD00[quadPt][j] * Mean_Comp2[globDOF];
                    U2x_Mean[quadPt] += origvaluesD10[quadPt][j] * Mean_Comp2[globDOF];
                    U2y_Mean[quadPt] += origvaluesD01[quadPt][j] * Mean_Comp2[globDOF];
                }
            }
            double U1_Mode_a[N_Points2];
            double U2_Mode_a[N_Points2];

            double U1x_Mode_a[N_Points2];
            double U1y_Mode_a[N_Points2];
            double U1xx_Mode_a[N_Points2];
            double U1yy_Mode_a[N_Points2];

            double U2x_Mode_a[N_Points2];
            double U2y_Mode_a[N_Points2];
            double U2xx_Mode_a[N_Points2];
            double U2yy_Mode_a[N_Points2];

            for (int quadPt = 0; quadPt < N_Points2; quadPt++) // Initialize
            {
                // C_i[quadPt] = 0;
                U1_Mode_a[quadPt] = 0;

                U2_Mode_a[quadPt] = 0;

                U1x_Mode_a[quadPt] = 0;
                U1y_Mode_a[quadPt] = 0;
                U1xx_Mode_a[quadPt] = 0;
                U1yy_Mode_a[quadPt] = 0;
                U2xx_Mode_a[quadPt] = 0;
                U2yy_Mode_a[quadPt] = 0;
            }

            for (int quadPt = 0; quadPt < N_Points2; quadPt++)
            {
                for (int j = 0; j < N_BaseFunct; j++)
                {
                    int globDOF = DOF[j];
                    U1_Mode_a[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp1_a[globDOF];

                    U2_Mode_a[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp2_a[globDOF];

                    U1x_Mode_a[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp1_a[globDOF];
                    U1y_Mode_a[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp1_a[globDOF];
                    U1xx_Mode_a[quadPt] += origvaluesD20[quadPt][j] * Mode_Comp1_a[globDOF];
                    U1yy_Mode_a[quadPt] += origvaluesD02[quadPt][j] * Mode_Comp1_a[globDOF];

                    U2x_Mode_a[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp2_a[globDOF];
                    U2y_Mode_a[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp2_a[globDOF];

                    U2xx_Mode_a[quadPt] += origvaluesD20[quadPt][j] * Mode_Comp2_a[globDOF];
                    U2yy_Mode_a[quadPt] += origvaluesD02[quadPt][j] * Mode_Comp2_a[globDOF];
                }
            }

            for (int qdpt = 0; qdpt < N_Points2; qdpt++)
            {
                double nu = Coeffs[qdpt][0];
                double Mult = Weights2[qdpt] * AbsDetjk[qdpt];
                ipval1 += (U1_Mode_a[qdpt] * U1x_Mean[qdpt] + U2_Mode_a[qdpt] * U1y_Mean[qdpt] + U1_Mean[qdpt] * U1x_Mode_a[qdpt] + U2_Mean[qdpt] * U1y_Mode_a[qdpt] - nu * U1xx_Mode_a[qdpt] - nu * U1yy_Mode_a[qdpt]) * U1_Mode_i[qdpt] * Mult;

                ipval1 += (U1_Mode_a[qdpt] * U2x_Mean[qdpt] + U2_Mode_a[qdpt] * U2y_Mean[qdpt] + U1_Mean[qdpt] * U2x_Mode_a[qdpt] + U2_Mean[qdpt] * U2y_Mode_a[qdpt] - nu * U2xx_Mode_a[qdpt] - nu * U2yy_Mode_a[qdpt]) * U2_Mode_i[qdpt] * Mult;
            }
            for (int i = 0; i < MaxN_QuadPoints_2D; i++)
            {
                delete[] Coeffs[i];
            }
            delete[] Coeffs;
        } // cell loop ip1 end

        for (int r = 0; r < lenPhi; r++)
        {
            phi_New[r] += (ipval1 * phi_Array_a[r] * -1.0);
        }
        ipval1 = 0.0;

    } // a loop end ip1
    for (int a = 0; a < N_S; a++)
    { // a loop ip2
        for (int b = 0; b < N_S; b++)
        { // b loop ip2
            for (int cellId = 0; cellId < N_Cells; cellId++)
            { // cell loop ip2
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
                    int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
                    QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
                    TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
                    QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

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
                    int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
                    QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
                    TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
                    QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

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
                val = 0;

                // Get Coefficients b1 and b2
                double *Param[MaxN_QuadPoints_2D];
                double **Coeffs = new double *[MaxN_QuadPoints_2D];
                for (int i = 0; i < MaxN_QuadPoints_2D; i++)
                {
                    Coeffs[i] = new double[10]();
                }
                DO_Mode_Equation_Coefficients(N_Points2, X, Y, Param, Coeffs);
                // Save Values of C at all quadrature points for I component
                memcpy(Mode_Comp1_a, U_Mode + (a * 2 * lenMode), lenMode * SizeOfDouble);       // col Major
                memcpy(Mode_Comp2_a, U_Mode + ((a * 2 + 1) * lenMode), lenMode * SizeOfDouble); // col Major
                memcpy(Mode_Comp1_b, U_Mode + (b * 2 * lenMode), lenMode * SizeOfDouble);       // col Major
                memcpy(Mode_Comp2_b, U_Mode + ((b * 2 + 1) * lenMode), lenMode * SizeOfDouble); // col Major
                memcpy(phi_Array_a, Phi_Array + a * lenPhi, lenPhi * SizeOfDouble);
                memcpy(phi_Array_b, Phi_Array + b * lenPhi, lenPhi * SizeOfDouble);

                double U1_Mode_i[N_Points2];
                double U2_Mode_i[N_Points2];
                double U1_Mode_a[N_Points2];
                double U2_Mode_a[N_Points2];
                double U1x_Mode_b[N_Points2];
                double U1y_Mode_b[N_Points2];
                double U2x_Mode_b[N_Points2];
                double U2y_Mode_b[N_Points2];

                for (int quadPt = 0; quadPt < N_Points2; quadPt++) // Initialize
                {
                    // C_i[quadPt] = 0;
                    U1_Mode_i[quadPt] = 0;
                    U2_Mode_i[quadPt] = 0;
                    U1_Mode_a[quadPt] = 0;
                    U2_Mode_a[quadPt] = 0;
                    U1x_Mode_b[quadPt] = 0;
                    U1y_Mode_b[quadPt] = 0;
                    U2x_Mode_b[quadPt] = 0;
                    U2y_Mode_b[quadPt] = 0;
                }

                for (int quadPt = 0; quadPt < N_Points2; quadPt++)
                {
                    for (int j = 0; j < N_BaseFunct; j++)
                    {
                        int globDOF = DOF[j];
                        U1_Mode_i[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp1_i[globDOF];
                        U2_Mode_i[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp2_i[globDOF];
                        U1_Mode_a[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp1_a[globDOF];
                        U2_Mode_a[quadPt] += origvaluesD00[quadPt][j] * Mode_Comp2_a[globDOF];
                        U1x_Mode_b[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp1_b[globDOF];
                        U1y_Mode_b[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp1_b[globDOF];
                        U2x_Mode_b[quadPt] += origvaluesD10[quadPt][j] * Mode_Comp2_b[globDOF];
                        U2y_Mode_b[quadPt] += origvaluesD01[quadPt][j] * Mode_Comp2_b[globDOF];
                    }
                }
                for (int qdpt = 0; qdpt < N_Points2; qdpt++)
                {
                    double Mult = Weights2[qdpt] * AbsDetjk[qdpt];

                    ipval2 += (U1_Mode_a[qdpt] * U1x_Mode_b[qdpt] + U2_Mode_a[qdpt] * U1y_Mode_b[qdpt]) * U1_Mode_i[qdpt] * Mult;
                    ipval2 += (U1_Mode_a[qdpt] * U2x_Mode_b[qdpt] + U2_Mode_a[qdpt] * U2y_Mode_b[qdpt]) * U2_Mode_i[qdpt] * Mult;
                }
                for (int i = 0; i < MaxN_QuadPoints_2D; i++)
                {
                    delete[] Coeffs[i];
                }
                delete[] Coeffs;
            } // cell loop end ip2

            for (int r = 0; r < lenPhi; r++)
            {
                phi_New[r] += (-1.0 * (phi_Array_a[r] * phi_Array_b[r] - TDatabase::ParamDB->COVARIANCE_MATRIX_DO[a * N_S + b]) * ipval2);
            }
            ipval2 = 0.0;

        } // b loop end ip2
    }     // a loop end ip2
    cout << "phi new" << Ddot(lenPhi, phi_New, phi_New) << endl;
    double timeStep = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    for (int i = 0; i < lenPhi; i++)
    {
        phi_Array_i[i] = phi_Old_i[i] + timeStep * phi_New[i];
    }

    delete[] phi_New;
    delete[] Mean_Comp1;
    delete[] Mean_Comp2;
    delete[] Phi_Old;
    delete[] Mode_Comp1_i;
    delete[] Mode_Comp2_i;
    delete[] phi_Old_i;
    delete[] Mode_Comp1_a;
    delete[] Mode_Comp2_a;
    delete[] phi_Array_a;
    delete[] Mode_Comp1_b;
    delete[] Mode_Comp2_b;
    delete[] phi_Array_b;
}

void DO_Mode_RHS_Aux_Param(double *in, double *out)
{
    out = in;
}
