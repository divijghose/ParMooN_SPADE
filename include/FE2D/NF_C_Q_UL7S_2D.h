/*
TNodalFunctional2D(NodalFunctional2D id,
         int n_allfunctionals, int n_edgefunctionals, 
         int n_pointsall, int n_pointSdge, 
         double *xi, double *eta, double *t, 
         DoubleFunctVect *evalall, 
         DoubleFunctVect *evaledge); 
*/

static double NF_C_Q_UL7S_2D_Xi[] = {
-1.0 ,-5.0/7.0 ,-3.0/7.0 ,-1.0/7.0 ,1.0/7.0 ,3.0/7.0 ,5.0/7.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,5.0/7.0 ,3.0/7.0 ,1.0/7.0 ,-1.0/7.0 ,-3.0/7.0 ,-5.0/7.0 ,-1.0 ,-1.0 ,-1.0 ,-1.0 ,-1.0 ,-1.0 ,-1.0 ,-5.0/7.0 ,-5.0/7.0 ,-5.0/7.0 ,-5.0/7.0 ,-5.0/7.0 ,-3.0/7.0 ,-3.0/7.0 ,-3.0/7.0 ,-3.0/7.0 ,-3.0/7.0 ,-1.0/7.0 ,-1.0/7.0 ,-1.0/7.0 ,-1.0/7.0 ,-1.0/7.0 ,1.0/7.0 ,1.0/7.0 ,1.0/7.0 ,1.0/7.0 ,1.0/7.0 ,3.0/7.0 ,3.0/7.0 ,3.0/7.0 ,3.0/7.0 ,3.0/7.0 ,5.0/7.0 ,-5.0/7.0
};

static double NF_C_Q_UL7S_2D_Eta[] = {
-1.0 ,-1.0 ,-1.0 ,-1.0 ,-1.0 ,-1.0 ,-1.0 ,-1.0 ,-5.0/7.0 ,-3.0/7.0 ,-1.0/7.0 ,1.0/7.0 ,3.0/7.0 ,5.0/7.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,5.0/7.0 ,3.0/7.0 ,1.0/7.0 ,-1.0/7.0 ,-3.0/7.0 ,-5.0/7.0 ,-5.0/7.0 ,-3.0/7.0 ,-1.0/7.0 ,1.0/7.0 ,3.0/7.0 ,-5.0/7.0 ,-3.0/7.0 ,-1.0/7.0 ,1.0/7.0 ,3.0/7.0 ,-5.0/7.0 ,-3.0/7.0 ,-1.0/7.0 ,1.0/7.0 ,3.0/7.0 ,-5.0/7.0 ,-3.0/7.0 ,-1.0/7.0 ,1.0/7.0 ,3.0/7.0 ,-5.0/7.0 ,-3.0/7.0 ,-1.0/7.0 ,1.0/7.0 ,3.0/7.0 ,-5.0/7.0 ,5.0/7.0
};

static double NF_C_Q_UL7S_2D_T[] = {
-1.0 ,-5.0/7.0 ,-3.0/7.0 ,-1.0/7.0 ,1.0/7.0 ,3.0/7.0 ,5.0/7.0 ,1.0 };

void NF_C_Q_UL7S_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  memcpy(Functionals, PointValues, 55*SizeOfDouble);
};

void NF_C_Q_UL7S_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
  memcpy(Functionals, PointValues, 8*SizeOfDouble);
};

TNodalFunctional2D *NF_C_Q_UL7S_2D_Obj = new TNodalFunctional2D
        (NF_C_Q_UL7S_2D, 55, 8, 55, 8, NF_C_Q_UL7S_2D_Xi, NF_C_Q_UL7S_2D_Eta,
         NF_C_Q_UL7S_2D_T, NF_C_Q_UL7S_2D_EvalAll, NF_C_Q_UL7S_2D_EvalEdge);
