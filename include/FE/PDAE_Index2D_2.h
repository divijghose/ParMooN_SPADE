// ======================================================================
// @(#)TimeConvDiff2D.h        1.6 04/13/00
//
// common declaration for time dependent convection diffusion problems
// ======================================================================

#ifndef __PDAE_INDEX2D_2__
#define __PDAE_INDEX2D_2__

// part for standard Galerkin
int N_Terms = 3;
MultiIndex2D Derivatives[3] = { D10, D01, D00 };
int SpacesNumbers[3] = { 0, 0, 0 };

// part for SDFEM
int N_Terms_SD = 5;
MultiIndex2D Derivatives_SD[5] = { D10, D01, D00, D20, D02 };
int SpacesNumbers_SD[5] = { 0, 0, 0, 0 };

// part for all
int N_Matrices = 4;
int RowSpace[4] = { 0, 0, 0, 0  };
int ColumnSpace[4] = { 0, 0, 0, 0 };
int N_Rhs = 2;
int RhsSpace[2] = { 0, 0 };

// part for mass matrix
int N_Matrices_Mass = 4;
int RowSpace_Mass[4] = { 0, 0, 0, 0 };
int ColumnSpace_Mass[4] = { 0, 0, 0, 0 };
int N_Rhs_Mass= 0;
int N_Terms_Mass = 1;
MultiIndex2D Derivatives_Mass[1] = { D00 };
int SpacesNumbers_Mass[1] = { 0 };

MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };

void TimeMassAssemble_PDAE2(double Mult, double *coeff, double hK,
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs);

void TimeBilinearAssemble_PDAE2(double Mult, double *coeff, double *param, double hK,
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs);

void TimeBilinearAssembleJ_PDAE2(double Mult, double *coeff, double *param, double hK,
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs);


void TimeBilinearAssembleC_PDAE2(double Mult, double *coeff, double *param, double hK,
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs);

void TimeRhsAssemble_PDAE2(double Mult, double *coeff, double *param, double hK,
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs);

#endif
