// ***********************************************************************
// P1 element, discontinous, 2D, triangle
// ***********************************************************************

// number of degrees of freedom
static int D_T_P1_2D_NDOF = 3;

// number of dofs on the closure of the joints
static int D_T_P1_2D_JointDOF = 0;

// which local dofs are on the joints
static int *D_T_P1_2D_J0 = NULL;
static int *D_T_P1_2D_J1 = NULL;
static int *D_T_P1_2D_J2 = NULL;

static int *D_T_P1_2D_J[3] = { D_T_P1_2D_J0, D_T_P1_2D_J1, D_T_P1_2D_J2 };

// number of inner dofs
static int D_T_P1_2D_NInner = 3;

// array containing the numbers for the inner dofs (here is no inner dof)
static int D_T_P1_2D_Inner[3] = { 0, 1, 2 };

static char D_T_P1_2D_String[] = "D_T_P1_2D";

TFEDesc2D *FE_D_T_P1_2D_Obj=new TFEDesc2D(D_T_P1_2D_String, D_T_P1_2D_NDOF,
                              D_T_P1_2D_JointDOF, D_T_P1_2D_J,
                              D_T_P1_2D_NInner, D_T_P1_2D_Inner);
