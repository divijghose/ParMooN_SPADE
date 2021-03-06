// ======================================================================
// instationary problem
// ======================================================================

/// ========================================================================
// example file
// ========================================================================

void ExampleFile()
{
  OutPut("Example: SinCos4.h" << endl); 
  TDatabase::ParamDB->INTERNAL_STEADY_STATE_MATRICES_OR_RHS = 0;
}
// exact solution
void Exact(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*t*cos(x*y*y);
  values[1] = -(t*t)*sin(x*y*y)*y*y;
  values[2] = -(t*t)*sin(x*y*y)*2*x*y;
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;
   
  switch(BdComp)
  {
    case 0: value = t*t;
            break;
    case 1: value = t*t*cos(Param*Param);
            break;
    case 2: value = t*t*cos(1-Param);
            break;
    case 3: value = t*t;
            break;
  }
}

// initial conditon
void InitialCondition(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*t*cos(x*y*y);
}


void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR;
  double a=2, b=-1, c=1;
  int i;
  double *coeff, *param;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = a;
    coeff[2] = b;
    coeff[3] = c;

    coeff[4] = 2*t*cos(x*y*y)
       - eps* t*t*(-cos(x*y*y)*(y*y*y*y+4*x*x*y*y) - sin(x*y*y)*2*x)  
       - a*t*t*sin(x*y*y)*y*y - b*t*t*sin(x*y*y)*2*x*y 
       + c*t*t*cos(x*y*y);
  }
}



















