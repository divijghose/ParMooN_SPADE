/*
TNodalFunctional2D(NodalFunctional2D id,
         int n_allfunctionals, int n_edgefunctionals, 
         int n_pointsall, int n_pointsedge, 
         double *xi, double *eta, double *t, 
         DoubleFunctVect *evalall, 
         DoubleFunctVect *evaledge); 
*/

static double NF_C_Q_Q4_2D_Xi[] = {
        -1.000000000000000e+00,
        -5.000000000000000e-01,
         0.000000000000000e+00,
         5.000000000000000e-01,
         1.000000000000000e+00,
        -1.000000000000000e+00,
        -5.000000000000000e-01,
         0.000000000000000e+00,
         5.000000000000000e-01,
         1.000000000000000e+00,
        -1.000000000000000e+00,
        -5.000000000000000e-01,
         0.000000000000000e+00,
         5.000000000000000e-01,
         1.000000000000000e+00,
        -1.000000000000000e+00,
        -5.000000000000000e-01,
         0.000000000000000e+00,
         5.000000000000000e-01,
         1.000000000000000e+00,
        -1.000000000000000e+00,
        -5.000000000000000e-01,
         0.000000000000000e+00,
         5.000000000000000e-01,
         1.000000000000000e+00
};

static double NF_C_Q_Q4_2D_Eta[] = {
        -1.000000000000000e+00,
        -1.000000000000000e+00,
        -1.000000000000000e+00,
        -1.000000000000000e+00,
        -1.000000000000000e+00,
        -5.000000000000000e-01,
        -5.000000000000000e-01,
        -5.000000000000000e-01,
        -5.000000000000000e-01,
        -5.000000000000000e-01,
         0.000000000000000e+00,
         0.000000000000000e+00,
         0.000000000000000e+00,
         0.000000000000000e+00,
         0.000000000000000e+00,
         5.000000000000000e-01,
         5.000000000000000e-01,
         5.000000000000000e-01,
         5.000000000000000e-01,
         5.000000000000000e-01,
         1.000000000000000e+00,
         1.000000000000000e+00,
         1.000000000000000e+00,
         1.000000000000000e+00,
         1.000000000000000e+00
};

static double NF_C_Q_Q4_2D_T[] = {
        -1.000000000000000e+00,
        -5.000000000000000e-01,
         0.000000000000000e+00,
         5.000000000000000e-01,
         1.000000000000000e+00
};

void NF_C_Q_Q4_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  memcpy(Functionals, PointValues, 25*SizeOfDouble);
};

void NF_C_Q_Q4_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
  memcpy(Functionals, PointValues, 5*SizeOfDouble);
};

TNodalFunctional2D *NF_C_Q_Q4_2D_Obj = new TNodalFunctional2D
        (NF_C_Q_Q4_2D, 25, 5, 25, 5, NF_C_Q_Q4_2D_Xi, NF_C_Q_Q4_2D_Eta,
         NF_C_Q_Q4_2D_T, NF_C_Q_Q4_2D_EvalAll, NF_C_Q_Q4_2D_EvalEdge);
