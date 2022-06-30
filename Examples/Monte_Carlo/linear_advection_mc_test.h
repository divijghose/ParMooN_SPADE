
/**
 * @file linear_advection_mc_test.h
 * @brief Purpose:     Example file for solving the set of dynamically orthogonal field
                       equations for linear advection.
                       Features included in this example file :
                       1. Definition of boundary conditions, boundary values,
                       initial condition and bilinear coefficients
                       2. Assembly functions for Mean Equation, Mode Equation
                       and Coefficient Equation
                       3. Assembly function for RHS of Mode Equation

 * @authors Sashikumaar Ganesan
 * @authors Divij Ghose
 * @authors Thivin Anandh

 * @bug No known bugs
 *
 */

// ===========================================================================//
//  Dynamically Orthogonal Field Equation Solution of Linear Advection Problem //
// ===========================================================================//

// =======================================================================
//
// Purpose:     Example file for solving the set of dynamically orthogonal field
//              equations for linear advection.
//              Features included in this example file -
//              1. Definition of boundary conditions, boundary values,
//                 initial condition and bilinear coefficients
//              2. Assembly functions for Mean Equation, Mode Equation
//                 and Coefficient Equation
//              3. Assembly function for RHS of Mode Equation
//
// Authors:      Sashikumaar Ganesan, Thivin Anandh, Divij Ghose
//
// History:     1> First iteration implemented on 18.03.2022
//				2> Bug fixes on 24.03.2022

// =======================================================================

#include <MacroCell.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>
#include <memory>

extern "C"
{
#include <gridgen.h>

    void triangulate(char *, struct triangulateio *,
                     struct triangulateio *, struct triangulateio *);
}

// Include files related to the cell looping for the RHS filling part.

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <FESpace2D.h>
#include <SystemTNSE2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <Output2D.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>
// #include <TimeUtilities.h>
#include <TNSE2D_ParamRout.h>
#include <TimeDiscRout.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <QuadAffin.h>
#include <QuadBilinear.h>

//========================================================================//
// 								Example File							  //
// =======================================================================//

#define __SIN3__

void ExampleFile()
{
    OutPut("Example: Linear Advection - Dynamically Orthogonal Field Equation Solution" << endl);
}

// exact solution

void Exact(double x, double y, double *values)
{
    double t;

    t = TDatabase::TimeDB->CURRENTTIME;

    values[0] = exp(t) * (sin(2 * Pi * x) * sin(2 * Pi * y));
    values[1] = exp(t) * 2 * Pi * cos(2 * Pi * x) * sin(2 * Pi * y);
    values[2] = exp(t) * 2 * Pi * sin(2 * Pi * x) * cos(2 * Pi * y);
    values[3] = 0;
}

// kind of boundary condition (for FE space needed)
/**
 * @brief Boundary Conditions for the linear advection problem. For a unit square geometry, using ParMooN's internal mesh generator, the boundary ID's are generated by default, starting from 0 at the bottom boundary and increasing in an anti-clockwise manner.
 * @param BdComp Boundary ID, integer number used to identify the boundary of the geometry
 * @param t Param value
 * @param cond Boundary Condition at given BdComp (i.e. DIRICHLET or NEUMANN)
 * @remark For the Linear Advection DO template problem, we apply a DIRICHLET boundary condition on all 4 boundaries of the unit square domain.
 * @remark This test example file also includes the Monte Carlo solutions. The boundary condition applies to all Monte Carlo realizations.
 */
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
    //  cond = NEUMANN;

    //  if(BdComp == 1 || BdComp == 3 || BdComp == 2 )
    //  {
    cond = DIRICHLET;
    //  }
}

/**
 * @brief Value of prescribed boundary conditions for the linear advection problem. For a unit square geometry, using ParMooN's internal mesh generator, the boundary ID's are generated by default, starting from 0 at the bottom boundary and increasing in an anti-clockwise manner.
 *
 * @param BdComp Boundary ID, integer number used to identify the boundary of the geometry
 * @param Param Param value
 * @param value Value of the prescribed boundary condition at given boundary
 *
 * @remark For the Linear Advection DO template problem, the value of the DIRICHLET boundary condition is set to 0 for all 4 boundaries of the unit square domain. This calue applied to both mean and the mode.
 * @remark This test example file also includes the Monte Carlo solutions. The boundary values applies to all Monte Carlo realizations.
 */
// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
    value = 0;
}

// initial conditon
/**
 * @brief Initial Condition for the linear advection problem.
 *
 * @param x X co-ordinate
 * @param y Y co-ordinate
 * @param values Value of initial condition
 */
void InitialCondition(double x, double y, double *values)
{
    double t, val;

    t = TDatabase::TimeDB->CURRENTTIME;
}
/**
 * @brief Bilinear Coefficients for the linear advection problem. Since the linear advection problem is solved on the template TCD2D problem, the diffusion, reaction and source parameters are zero, i.e.
 * coeff[0]=0 (Diffusion Parameter)
 * coeff[3]=0 (Reaction Parameter)
 * coeff[4]=0 (Source term)
 *
 *
 * @param n_points
 * @param X X co-ordinate
 * @param Y Y co-oridnate
 * @param parameters
 * @param coeffs Bilinear Coefficients
 *
 * @remark The X and Y components of the advection velocity, coeff[1] and coeff[2] are adjusted to give a circular, anti-clockwise advection field -
 * coeff[1] = -sin(t)*constant
 * coeff[2] = cos(t)*constant
 * @remark This function can be used to define BilinearCoeffs for any general problem.
 */
void BilinearCoeffs(int n_points, double *X, double *Y,
                    double **parameters, double **coeffs)
{
    double eps = 1 / TDatabase::ParamDB->PE_NR;
    double b1 = 1, b2 = -1, c = 1;
    int i;
    double *coeff;
    double x, y;
    double t = TDatabase::TimeDB->CURRENTTIME;

    for (i = 0; i < n_points; i++)
    {
        coeff = coeffs[i];

        x = X[i];
        y = Y[i];

        coeff[0] = 0;
        coeff[1] = -sin(t) * 0.2;
        coeff[2] = cos(t) * 0.2;
        coeff[3] = 0;

        coeff[4] = 0.0;
    }
}

void readRealizationFromText(double *RealznVect, const int N_R, const int N_DOF)
{
    cout << "Read In" << endl;
    std::vector<std::vector<std::string>> content;
    std::vector<std::string> row;
    std::string line, word;

    std::string fileInName = "Realizations_" + std::to_string(N_R) + "_NDOF_" + std::to_string(N_DOF) + ".txt";
    std::ifstream file(fileInName);
    if (file.is_open())
    {
        while (getline(file, line))
        {
            row.clear();

            std::stringstream str(line);

            while (getline(str, word, ','))
                row.push_back(word);
            content.push_back(row);
        }
        cout << "Realization file opened succesfully" << endl;
    }
    else
    {
        cout << "Could not open the file\n"
             << endl
             << "Please check is file " << fileInName << " exists" << endl;
        exit(0);
    }

    cout << "" << endl;
    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = 0; j < N_R; j++)
        {
            RealznVect[i * N_R + j] = std::stod(content[i][j]);
        }
    }

    cout << "Realization file read successfully" << endl;
    return;
}

void writeRealizationToText(const double *RealznVect, const int N_R, const int N_DOF)
{
    std::string fileoutMC = "Realizations_" + std::to_string(N_R) + "_NDOF_" + std::to_string(N_DOF) + ".txt";
    std::ofstream fileMC;
    fileMC.open(fileoutMC);

    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = 0; j < N_R; j++)
        {
            fileMC << RealznVect[i * N_R + j];
            if (j != N_R - 1)
                fileMC << ",";
        }
        fileMC << endl;
    }
    cout << "All Realizations written to: " << fileoutMC << endl;
    return;
}

void calcMeanRealization(const double *RealznVect, double *MeanVect, const int N_R, const int N_DOF)
{
    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = 0; j < N_R; j++)
        {
            MeanVect[i] += RealznVect[i * N_R + j] / N_R;
        }
    }
    cout << "Mean of realizations computed successfully" << endl;
    return;
}

void calcStdDevRealization(const double *RealznVect, double *StdDevVect, const int N_R, const int N_DOF)
{
    double *MeanVector = new double[N_DOF]();
    calcMeanRealization(RealznVect, MeanVector, N_R, N_DOF);
    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = 0; j < N_R; j++)
        {
            StdDevVect[i] += pow(RealznVect[i * N_R + j] - MeanVector[i], 2) / N_R;
        }

        StdDevVect[i] = sqrt(StdDevVect[i]);
    }

    delete[] MeanVector;
    cout << "Standard deviation of realizations computed successfully" << endl;


    return;
}