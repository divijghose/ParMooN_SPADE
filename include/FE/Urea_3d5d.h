/****************************************************************************************
 *                                                                                      *
 *                         Urea_3d5d.h                                                  *
 *                        -------------                                                 *
 *                                                                                      *
 ***************************************************************************************/

#ifndef __UREA_3D5D__
#define __UREA_3D5D__
double Volume_2D(int index,int  N_x, int  N_y,int  N_z, int  N_a, int  N_b, double *a_layers_coord, double *b_layers_coord,
                                         double *f);
double growth_rate_5D(double c, double temp, int variable);
double b_nuc_5D(double c, double temp);
double saturation_5D(double c, double temp);

double psd_boundary_urea_5D(double x, double y, double z, double a, double b,double *a_layers_coord,double *b_layers_coord,
                    double conc, double temp);
int PSD_bound_cond_from_velo_inflow_urea_5D(double x, double y, double z);
void PrepareNucleation(TFESpace3D *fespace,
TFEFunction3D *temperature,TFEFunction3D *concentration,
int N_x, int N_y, int N_z, int N_a,int N_b,
double *x_coord, double *y_coord, double *z_coord,
double *a_layers_coord, double *b_layers_coord,double *f, double *rhs);
void Urea_RKV_FDM_5D(TCollection *coll,
TFEFunction3D *velocity1, TFEFunction3D *velocity2, TFEFunction3D *velocity3,
TFEFunction3D *concent_C,
TFEFunction3D *Temp,
double *f_old, double *rhs_psd,  double **stages,
int N_x, int N_y, int N_z, int N_a,int N_b,
double *x_coord, double *y_coord, double *z_coord, double *a_coord, double *b_coord, double *a_layers_coord, double *b_layers_coord,
double *velo1, double *velo2, double *velo3, double *concent_C_array, double *Temp_array,
int *correspond_3dgrid);

void Calculate_Volume_Distribution_5D(TCollection *Coll,int N_x, int N_y, int N_z, int N_a, 
			   double *x_coord, double *y_coord, double *z_coord,
				   double *a_layers_coord, double *f);

void Calculate_PSD_on_node_5D(int N_x, int N_y, int N_z, int N_a, int N_b,
double *x_coord, double *y_coord, double *z_coord,
			   double *a_layers_coord, double *b_layers_coord, double *sol_psd, double x, double y, double z);


double psd_value_5D(double *f, double *a_layers_coord, double val, int index, int N_a, int N3);

#endif


