// ***********************************************************************
// P2 element, discontinuous, 3D Tetrahedra
// 
// Author:     Markus Wolff
//
// ***********************************************************************

static double D_T_P2_3D_CM[100] = {
    600,-1800,-1800,-1800,1260,2520,2520,1260,2520,1260,
    -1800,8640,4320,4320,-7560,-10080,-10080,-2520,-5040,-2520,
    -1800,4320,8640,4320,-2520,-10080,-5040,-7560,-10080,-2520,
    -1800,4320,4320,8640,-2520,-5040,-10080,-2520,-10080,-7560,
    1260,-7560,-2520,-2520,7560,7560,7560,1260,2520,1260,
    2520,-10080,-10080,-5040,7560,20160,10080,7560,10080,2520,
    2520,-10080,-5040,-10080,7560,10080,20160,2520,10080,7560,
    1260,-2520,-7560,-2520,1260,7560,2520,7560,7560,1260,
    2520,-5040,-10080,-10080,2520,10080,10080,7560,20160,7560,
    1260,-2520,-2520,-7560,1260,2520,7560,1260,7560,7560
};

static void D_T_P2_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  int nBF = 10; // number of basis functions
  double mon[]={1,xi,eta,zeta,xi*xi,xi*eta,xi*zeta,
                eta*eta,eta*zeta,zeta*zeta};
  memset(values, 0.0, nBF*SizeOfDouble);
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i] += D_T_P2_3D_CM[i+j*nBF]*mon[j];

    }
  }
}

static void D_T_P2_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  int nBF = 10; // number of basis functions
  double mon[]={0,1,0,0,2*xi,eta,zeta,
                0,0,0};
  memset(values, 0.0, nBF*SizeOfDouble);
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i] += D_T_P2_3D_CM[i+j*nBF]*mon[j];

    }
  }
}

static void D_T_P2_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
  int nBF = 10; // number of basis functions
  double mon[]={0,0,1,0,0,xi,0,
                2*eta,zeta,0};
  memset(values, 0.0, nBF*SizeOfDouble);
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i] += D_T_P2_3D_CM[i+j*nBF]*mon[j];

    }
  }
}

static void D_T_P2_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  int nBF = 10; // number of basis functions
  double mon[]={0,0,0,1,0,0,xi,
                0,eta,2*zeta};
  memset(values, 0.0, nBF*SizeOfDouble);
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i] += D_T_P2_3D_CM[i+j*nBF]*mon[j];

    }
  }
}

static void D_T_P2_3D_DeriveXiXi(double xi, double eta, double zeta,
                             double *values)
{
  int nBF = 10; // number of basis functions
  double mon[]={0,0,0,0,2,0,0,
                0,0,0};
  memset(values, 0.0, nBF*SizeOfDouble);
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i] += D_T_P2_3D_CM[i+j*nBF]*mon[j];

    }
  }
}

static void D_T_P2_3D_DeriveXiEta(double xi, double eta, double zeta,
                             double *values)
{
  int nBF = 10; // number of basis functions
  double mon[]={0,0,0,0,0,1,0,
                0,0,0};
  memset(values, 0.0, nBF*SizeOfDouble);
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i] += D_T_P2_3D_CM[i+j*nBF]*mon[j];

    }
  }
}

static void D_T_P2_3D_DeriveXiZeta(double xi, double eta, double zeta,
                             double *values)
{
  int nBF = 10; // number of basis functions
  double mon[]={0,0,0,0,0,0,1,
                0,0,0};
  memset(values, 0.0, nBF*SizeOfDouble);
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i] += D_T_P2_3D_CM[i+j*nBF]*mon[j];

    }
  }
}

static void D_T_P2_3D_DeriveEtaEta(double xi, double eta, double zeta,
                             double *values)
{
  int nBF = 10; // number of basis functions
  double mon[]={0,0,0,0,0,0,0,
                2,0,0};
  memset(values, 0.0, nBF*SizeOfDouble);
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i] += D_T_P2_3D_CM[i+j*nBF]*mon[j];

    }
  }
}

static void D_T_P2_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                             double *values)
{
  int nBF = 10; // number of basis functions
  double mon[]={0,0,0,0,0,0,0,
                0,1,0};
  memset(values, 0.0, nBF*SizeOfDouble);
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i] += D_T_P2_3D_CM[i+j*nBF]*mon[j];

    }
  }
}

static void D_T_P2_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                             double *values)
{
  int nBF = 10; // number of basis functions
  double mon[]={0,0,0,0,0,0,0,
                0,0,2};
  memset(values, 0.0, nBF*SizeOfDouble);
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i] += D_T_P2_3D_CM[i+j*nBF]*mon[j];

    }
  }
}

TBaseFunct3D *BF_D_T_P2_3D_Obj =
new TBaseFunct3D(10, BF_D_T_P2_3D, BFUnitTetrahedron,
                 D_T_P2_3D_Funct, D_T_P2_3D_DeriveXi,
                 D_T_P2_3D_DeriveEta, D_T_P2_3D_DeriveZeta,
                 D_T_P2_3D_DeriveXiXi, D_T_P2_3D_DeriveXiEta,
                 D_T_P2_3D_DeriveXiZeta, D_T_P2_3D_DeriveEtaEta,
                 D_T_P2_3D_DeriveEtaZeta, D_T_P2_3D_DeriveZetaZeta,
                 2, 1,
		 0, NULL);
