/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

static double NF_N_Q_Q1_2D_Xi[] = 
        { -0.77459666924148337703585307995647992, 0, 
           0.77459666924148337703585307995647992, 
           1, 1, 1, 
           0.77459666924148337703585307995647992, 0, 
          -0.77459666924148337703585307995647992, 
          -1, -1, -1 };
static double NF_N_Q_Q1_2D_Eta[] = 
        { -1, -1, -1, 
          -0.77459666924148337703585307995647992, 0,
           0.77459666924148337703585307995647992,
           1, 1, 1,
           0.77459666924148337703585307995647992, 0,
          -0.77459666924148337703585307995647992 };
static double NF_N_Q_Q1_2D_T[] = 
        { -0.77459666924148337703585307995647992, 0,
           0.77459666924148337703585307995647992 };

void NF_N_Q_Q1_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  static double weights[3] = { 0.5555555555555555555555555555555556,
                               0.88888888888888888888888888888888889,
                               0.5555555555555555555555555555555556 };
  Functionals[0] = ( weights[0]*PointValues[0]
                    +weights[1]*PointValues[1]
                    +weights[2]*PointValues[2]) * 0.5;
  Functionals[1] = ( weights[0]*PointValues[3]
                    +weights[1]*PointValues[4]
                    +weights[2]*PointValues[5]) * 0.5;
  Functionals[2] = ( weights[0]*PointValues[6]
                    +weights[1]*PointValues[7]
                    +weights[2]*PointValues[8]) * 0.5;
  Functionals[3] = ( weights[0]*PointValues[9]
                    +weights[1]*PointValues[10]
                    +weights[2]*PointValues[11]) * 0.5;
}

void NF_N_Q_Q1_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
  static double weights[3] = { 0.5555555555555555555555555555555556,
                               0.88888888888888888888888888888888889,
                               0.5555555555555555555555555555555556 };
  Functionals[0] =(  weights[0]*PointValues[0]
                   +weights[1]*PointValues[1]
                   +weights[2]*PointValues[2])*0.5;
}

TNodalFunctional2D *NF_N_Q_Q1_2D_Obj = new TNodalFunctional2D
        (NF_N_Q_Q1_2D, 4, 1, 12, 3, NF_N_Q_Q1_2D_Xi, NF_N_Q_Q1_2D_Eta,
         NF_N_Q_Q1_2D_T, NF_N_Q_Q1_2D_EvalAll, NF_N_Q_Q1_2D_EvalEdge);
