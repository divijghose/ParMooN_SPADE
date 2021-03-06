// nonconforming P1 and Q1 like behaviour on considered edge
static char N5_2_N4_2_1Reg_Name[] = "N5_2_N5_2_1Reg";
static char N5_2_N4_2_1Reg_Desc[] = "nonconforming P5 or Q5 element, one regular grid";
static int N5_2_N4_2_1Reg_N0 = 5;
static int N5_2_N4_2_1Reg_N1 = 4;
static int N5_2_N4_2_1Reg_N2 = 4;
static int N5_2_N4_2_1Reg_NMid = 0;
static int *N5_2_N4_2_1Reg_Mid = NULL;
static int N5_2_N4_2_1Reg_NPairs = 0;
static int *N5_2_N4_2_1Reg_Pairs = NULL;
static int N5_2_N4_2_1Reg_NHanging = 4;
static int N5_2_N4_2_1Reg_Hanging[] = { 0, 1, 2, 3 };
static HNDesc N5_2_N4_2_1Reg_HangingTypes[] = { HN_N_P1_2D_0, HN_N_P2_2D_0,
                                                HN_N_P3_2D_0, HN_N_P4_2D_0 };
static int N5_2_N4_2_1Reg_Coupling_0[] = { 5,          9 };
static int N5_2_N4_2_1Reg_Coupling_1[] = { 5, 6,       9, 10 };
static int N5_2_N4_2_1Reg_Coupling_2[] = { 5, 6, 7,    9, 10, 11 };
static int N5_2_N4_2_1Reg_Coupling_3[] = { 5, 6, 7, 8, 9, 10, 11, 12 };
static int *N5_2_N4_2_1Reg_Coupling[] = { N5_2_N4_2_1Reg_Coupling_0,
                                          N5_2_N4_2_1Reg_Coupling_1,
                                          N5_2_N4_2_1Reg_Coupling_2,
                                          N5_2_N4_2_1Reg_Coupling_3 };
static int N5_2_N4_2_1Reg_NFarHanging = 0;
static int *N5_2_N4_2_1Reg_FarHanging = NULL;
static HNDesc *N5_2_N4_2_1Reg_FarHangingTypes = NULL;
static int ****N5_2_N4_2_1Reg_FarCoupling = NULL;
static int N5_2_N4_2_1Reg_NNoopposite = 9;
static int N5_2_N4_2_1Reg_Nopposite[] = { 4, 5, 6, 7, 8, 9, 10, 11, 12 };
static int N5_2_N4_2_1Reg_NNodes = 13;

TFE2DMapper1Reg *N5_2_N4_2_1Reg = new TFE2DMapper1Reg(
                N5_2_N4_2_1Reg_Name, N5_2_N4_2_1Reg_Desc,
                N5_2_N4_2_1Reg_N0, N5_2_N4_2_1Reg_N1, N5_2_N4_2_1Reg_N2,
                N5_2_N4_2_1Reg_NPairs, (int *)N5_2_N4_2_1Reg_Pairs,
                N5_2_N4_2_1Reg_NMid, (int *)N5_2_N4_2_1Reg_Mid,
                N5_2_N4_2_1Reg_NHanging, N5_2_N4_2_1Reg_Hanging,
                N5_2_N4_2_1Reg_HangingTypes, N5_2_N4_2_1Reg_Coupling,
                N5_2_N4_2_1Reg_NFarHanging, N5_2_N4_2_1Reg_FarHanging,
                N5_2_N4_2_1Reg_FarHangingTypes, N5_2_N4_2_1Reg_FarCoupling,
                N5_2_N4_2_1Reg_NNoopposite, N5_2_N4_2_1Reg_Nopposite,
                N5_2_N4_2_1Reg_NNodes);
