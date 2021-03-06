/** ************************************************************************ 
*
* @class     TSystemNSE3D
* @brief     stores the information of a 3D NSE system matrix 
* @author    Sashikumaar Ganesan, 
* @date      27.01.15
* @History    
 ************************************************************************  */


#ifndef __SYSTEMNSE3D__
#define __SYSTEMNSE3D__

#include <SquareMatrix3D.h>
#include <FEVectFunct3D.h>
#include <NSE_MultiGrid.h>
#include <ItMethod.h>
#include <AssembleMat3D.h>

#ifdef _MPI
//#include "mpi.h"
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#include <ParDirectSolver.h>
#endif

#ifdef _OMPONLY
#include <ParDirectSolver.h>
#endif

/**class for 3D  NSE system matrix */
class TSystemNSE3D
{
  protected:

#ifdef _MPI
    TParFEMapper3D **ParMapper_U, **ParMapper_P;
    TParFECommunicator3D **ParComm_U, **ParComm_P;
    MPI_Comm Comm;
    TParDirectSolver *DS;
#endif
    
#ifdef _OMPONLY
    TParDirectSolver *DS;
#endif
    
    /** Number of multigrid levels */
    int N_Levels;
    
    /** starting level for the solver, e.g. Start_Level = N_Levels-1 for direct solver */
    int Start_Level;  
    
    /** DOFs of velocity and pressure spaces */    
    int N_TotalDOF, N_U, N_P, N_Active, N_DirichletDof;
    
    /** Number of free surface/interface faces */
    int *N_FreeSurfFaces;
    
    /** Cell numbers and Joint numbers of free surface/interface faces */    
    int **FreeSurfCellNumbers, **FreeSurfJointNumbers;
    
    /** velocity and pressure fespace */
    TFESpace3D **U_Space, **P_Space, *FeSpaces[5];
    
    /** velo FE function */
    TFEVectFunct3D **Velocity;
    
    /** Fe functions of NSE */
    TFEFunction3D **Pressure;    

    double **SolArray, **RhsArray, *RHSs[4];

    TFESpace3D *fesp[4], *fesprhs[3], *fesp_aux[3];
    TFEFunction3D *fefct[7], **fefct_aux;
  
    // Assembling */
    TAssembleMat3D **AMatRhsAssemble, **AMatAssembleNonLinear;
    
    /** Discretization type */
    int Disctype;
       
    /** NSE type */
    int NSEType; 
           
    /** Bilinear coefficient   */
    CoeffFct3D *LinCoeffs[1];    
    
    /** NSEaux is used to pass additional fe functions (eg. mesh velocity) that is nedded for assembling */
    TAuxParam3D **NSEaux, **NSEaux_error;
    
    /** method for matrix vector mult */
    MatVecProc *MatVect;  
    
    /** method for resudual calculation */
    DefectProc *Defect;   
    
    /** Solver type */
    int SOLVER;
       
    /** number of matrices in the system matrix*/
    int N_Matrices;

    /** sqstructureA of the system matrix */
    TSquareStructure3D **sqstructureA;

    /** structure of the system matrix */
    TStructure3D **structureB, **structureBT;
    
    /** A is the stiffness/system mat for NSE velocity component   */
    TSquareMatrix3D **SqmatrixA11, **SqmatrixA12, **SqmatrixA13;
    TSquareMatrix3D **SqmatrixA21, **SqmatrixA22, **SqmatrixA23;
    TSquareMatrix3D **SqmatrixA31, **SqmatrixA32, **SqmatrixA33, *SQMATRICES[15];
    TSquareMatrix **sqmatrices;
    TSquareMatrix3D **SqmatrixF11, **SqmatrixF22, **SqmatrixF33;
    
    /** B is the  system mat for NSE pressure component   */
    TMatrix3D **MatrixB1, **MatrixB2, **MatrixB3, **MatrixB1T, **MatrixB2T, **MatrixB3T, *MATRICES[12];
    TMatrix **matrices;
    
    /** Boundary conditon */
    BoundCondFunct3D *BoundaryConditions[3];
  
     /** Boundary values */   
    BoundValueFunct3D *BoundaryValues[3];
        
    /** Discrete form for the equation */
    TDiscreteForm3D *DiscreteFormARhs, *DiscreteFormNL;

    /** variables for multigrid */
    int N_aux;
    double Parameters[2], *Itmethod_sol, *Itmethod_rhs;
    TNSE_MultiGrid *MG;
    TNSE_MGLevel *MGLevel;
    TItMethod *Itmethod, *prec;
 
    
  private:
   void UpdateUpwind(int i);
    
//    void UpdateLPS();
   
  public:
    /** constructor */
     TSystemNSE3D(int N_levels, TFESpace3D **velocity_fespace, TFESpace3D **presssure_fespace, TFEVectFunct3D **velocity, 
                     TFEFunction3D **pressure, double **sol, double **rhs, int disctype, int nsetype, int solver);

    /** destrcutor */
//     ~TSystemMatNSE3D();

    /** methods */
    /** Initilize the discrete forms and the matrices */    
    void Init(CoeffFct3D *lincoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *U1BoundValue, 
              BoundValueFunct3D *U2BoundValue, BoundValueFunct3D *U3BoundValue);
    
    /** assemble the system matrix */
    void Assemble();
    
    /** assemble the nonlinear part of the NSE system */
    void AssembleNonLinear(double **sol, double **rhs);
    
    /** get the resudual of the NSE system */
    void GetResidual(double *sol, double *rhs, double *res, double &impuls_residual, double &residual);
    
    /** solve the system matrix */
    void  Solve(double *sol, double *rhs);
  
    /** measure the error in the NSE */
    void MeasureErrors(DoubleFunct3D *ExactU1, DoubleFunct3D *ExactU2, DoubleFunct3D *ExactU3, DoubleFunct3D *ExactP,
                       double *u_error, double *p_error);
    
    /** find all joints which form the free surface/interface */   
    void FindFreeSurfJoints(int level, int Region_ID);
};

#endif
