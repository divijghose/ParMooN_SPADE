/**
 * @file burgers_mc.h
 * @brief Purpose:     Example file for Monte Carlo runs of Burgers' equation.
					   Features included in this example file :
					   1. Definition of boundary conditions, boundary values,
					   initial condition and bilinear coefficients
					   

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
// Purpose:     Example file for Monte Carlo runs of Burgers' equation.
//              Features included in this example file -
//              1. Definition of boundary conditions, boundary values,
//                 initial condition and bilinear coefficients
//      
//
// Authors:      Sashikumaar Ganesan, Thivin Anandh, Divij Ghose
//
// History:     1> First iteration implemented on 01.04.2022
//		        

// =======================================================================

void ExampleFile()
{

	OutPut("Example: Monte Carlo Solution for Time Dependent 2D Burgers' Equation" << endl);
}

// ========================================================================
// initial solution
// ========================================================================

void InitialU1(double x, double y, double *values)
{
	double t = TDatabase::TimeDB->CURRENTTIME;

	values[0] = 0;
}

void InitialU2(double x, double y, double *values)
{
	double t = TDatabase::TimeDB->CURRENTTIME;
	values[0] = 0;

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

}

void U2BoundValue(int BdComp, double Param, double &value)
{
	double t = TDatabase::TimeDB->CURRENTTIME, t1 = exp(-t);
	value = 0;
	
}

// data on each quadrature point for mean
void LinCoeffs(int n_points, double *X, double *Y,
			   double **parameters, double **coeffs)
{
	double nu = 1.0/TDatabase::ParamDB->RE_NR;
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
        coeff[0] = nu; // D(u):D(v) term
		coeff[1] = 0;
		coeff[2] = 0;
	}
}

