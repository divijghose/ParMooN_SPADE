// Navier-Stokes problem, exact solution is
// 
// u(x,y) = (4*y*(1-y),0)
// p(x,y) = 8*nu*(1/2-x)

void ExampleFile()
{
  OutPut("Example: PoiseuilleNu.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = 4*y*(1-y);
  values[1] = 0;
  values[2] = 4-8*y;
  values[3] = -8;
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
  double nu=1/TDatabase::ParamDB->RE_NR;
  values[0] = nu*(4-8*x);
  values[1] = -8*nu;
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
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=4*Param*(1-Param);
            break;
    case 2: value=0;
            break;
    case 3: value=4*Param*(1-Param);
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}


// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
  }
}


