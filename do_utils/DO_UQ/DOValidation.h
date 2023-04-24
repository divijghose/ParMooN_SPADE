#include <Database.h>
#include <FESpace2D.h>
#include <FEVectFunct2D.h>
#include <SystemTCD2D.h>
#include <Output2D.h>
#include <TimeDiscRout.h>

#include <mkl.h>
#include "Stats.h"

void validateDOvsMC(TFESpace2D *Scalar_FeSpace, double *RealizationVector, TOutput2D *Output, TSystemTCD2D *SystemMatrix)
{
    int N_Realisations = TDatabase::ParamDB->REALIZATIONS;
    int N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
    int subDim = TDatabase::ParamDB->N_Subspace_Dim;

    std::ostringstream os;
    os << " ";

    double *RealizationVectorCopy = new double[N_DOF * N_Realisations]();
    memcpy(RealizationVectorCopy, RealizationVector, N_DOF * N_Realisations * SizeOfDouble);

    double *MeanVectorMC = new double[N_DOF * 1]();
    calcMeanRealization(RealizationVectorCopy, MeanVectorMC, N_Realisations, N_DOF);

    double *stdDevVectorMC = new double[N_DOF * 1]();
    calcStdDevRealization(RealizationVectorCopy, stdDevVectorMC, N_Realisations, N_DOF);

    double *CompositeVectorMC = new double[N_DOF * 7]();
    for (int i = 0; i < N_DOF; i++)
    {
        CompositeVectorMC[i] = MeanVectorMC[i];
        CompositeVectorMC[N_DOF + i] = MeanVectorMC[i] + stdDevVectorMC[i];
        CompositeVectorMC[2 * N_DOF + i] = MeanVectorMC[i] - stdDevVectorMC[i];
        CompositeVectorMC[3 * N_DOF + i] = MeanVectorMC[i] + (2 * stdDevVectorMC[i]);
        CompositeVectorMC[4 * N_DOF + i] = MeanVectorMC[i] - (2 * stdDevVectorMC[i]);
        CompositeVectorMC[5 * N_DOF + i] = MeanVectorMC[i] + (3 * stdDevVectorMC[i]);
        CompositeVectorMC[6 * N_DOF + i] = MeanVectorMC[i] - (3 * stdDevVectorMC[i]);
    }

    delete[] MeanVectorMC;
    delete[] stdDevVectorMC;

    double *sol = new double[N_DOF]();
    double *rhs = new double[N_DOF]();
    double *oldrhs = new double[N_DOF]();
    TFEFunction2D *Scalar_FeFunction = new TFEFunction2D(Scalar_FeSpace, (char *)"sol", (char *)"sol", sol, N_DOF);

    Output->AddFEFunction(Scalar_FeFunction);

    int N_Composite = 7;
    int *imgVtk = new int[7]();
    int *imgTxt = new int[7]();
    std::string fileoutMC;

    for (int RealNo = 0; RealNo < N_Composite; RealNo++)
    {
        std::string filename;
        if (TDatabase::ParamDB->DOVerbose==1)
        {
            cout << " ============================================================================================================= " << endl;
            switch (RealNo)
            {
            case 0:
                cout << "Solving for Mean Solution" << endl;
                filename = "MonteCarlo_Mean_NR_" + std::to_string(N_Realisations);
                break;
            case 1:
                cout << "Solving for Mean + sigma Solution" << endl;
                filename = "MonteCarlo_MeanPlusSigma_NR_" + std::to_string(N_Realisations);
                break;
            case 2:
                cout << "Solving for Mean - sigma Solution" << endl;
                filename = "MonteCarlo_MeanMinusSigma_NR_" + std::to_string(N_Realisations);
                break;
            case 3:
                cout << "Solving for Mean + 2*sigma Solution" << endl;
                filename = "MonteCarlo_MeanPlus2Sigma_NR_" + std::to_string(N_Realisations);
                break;
            case 4:
                cout << "Solving for Mean - 2*sigma Solution" << endl;
                filename = "MonteCarlo_MeanMinus2Sigma_NR_" + std::to_string(N_Realisations);
                break;
            case 5:
                cout << "Solving for Mean + 3*sigma Solution" << endl;
                filename = "MonteCarlo_MeanPlus3Sigma_NR_" + std::to_string(N_Realisations);
                break;
            case 6:
                cout << "Solving for Mean - 3*sigma Solution" << endl;
                filename = "MonteCarlo_MeanMinus3Sigma_NR_" + std::to_string(N_Realisations);
                break;
            }
            cout << " ============================================================================================================= " << endl;
        }

        // Scalar_FeFunction->Interpolate(InitialCondition);
        std::string str = std::to_string(RealNo);

        char *VtkBaseName = const_cast<char *>(filename.c_str());

        for (int i = 0; i < N_DOF; i++)
            sol[i] = CompositeVectorMC[i + RealNo * N_DOF];

        printVTKOutput(VtkBaseName, imgVtk + RealNo, Output);

        fileoutMC = generateFileName("MonteCarlo/" + filename, imgTxt[RealNo], N_Realisations);
        printToTxt(fileoutMC, sol, N_DOF, 1, 'C');
        imgTxt[RealNo]++;
        // assemble the system matrix with given aux, sol and rhs
        // aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
        // otherwise, just pass with NULL
        // TSystemTCD2D *SystemMatrix = new TSystemTCD2D(Scalar_FeSpace, GALERKIN, DIRECT);

        SystemMatrix->AssembleMRhs(NULL, sol, rhs);

        //======================================================================
        // time disc loop
        //======================================================================
        // parameters for time stepping scheme
        int m = 0;
        int N_SubSteps = GetN_SubSteps();
        double end_time = TDatabase::TimeDB->ENDTIME;

        bool UpdateStiffnessMat = TRUE; // check BilinearCoeffs in example file
        bool UpdateRhs = TRUE;          // check BilinearCoeffs in example file
        bool ConvectionFirstTime = TRUE;

        // time loop starts
        while (TDatabase::TimeDB->CURRENTTIME < end_time)
        {
            m++;
            TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

            for (int l = 0; l < N_SubSteps; l++) // sub steps of fractional step theta
            {
                SetTimeDiscParameters(1);

                if (m == 1)
                {
                    // OutPut("Theta1: " << TDatabase::TimeDB->THETA1 << endl);
                    // OutPut("Theta2: " << TDatabase::TimeDB->THETA2 << endl);
                    // OutPut("Theta3: " << TDatabase::TimeDB->THETA3 << endl);
                    // OutPut("Theta4: " << TDatabase::TimeDB->THETA4 << endl);
                }

                double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
                TDatabase::TimeDB->CURRENTTIME += tau;

                // OutPut(endl<< "CURRENT TIME: ");
                // OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

                // copy rhs to oldrhs
                memcpy(oldrhs, rhs, N_DOF * SizeOfDouble);

                // unless the stiffness matrix or rhs change in time, it is enough to
                // assemble only once at the begning
                if (UpdateStiffnessMat || UpdateRhs || ConvectionFirstTime)
                {
                    SystemMatrix->AssembleARhs(NULL, sol, rhs);

                    // M:= M + (tau*THETA1)*A
                    // rhs: =(tau*THETA4)*rhs +(tau*THETA3)*oldrhs +[M-(tau*THETA2)A]*oldsol
                    // note! sol contains only the previous time step value, so just pass
                    // sol for oldsol
                    SystemMatrix->AssembleSystMat(oldrhs, sol, rhs, sol);
                    ConvectionFirstTime = FALSE;
                }

                // solve the system matrix
                SystemMatrix->Solve(sol, rhs);

                // restore the mass matrix for the next time step
                // unless the stiffness matrix or rhs change in time, it is not necessary to assemble the system matrix in every time step
                if (UpdateStiffnessMat || UpdateRhs)
                {
                    SystemMatrix->RestoreMassMat();
                }
                printVTKOutput(VtkBaseName, imgVtk + RealNo, Output);

                // if (TDatabase::ParamDB->WRITE_VTK)
                // {
                // os.seekp(std::ios::beg);
                // if (imgm[RealNo] < 10)
                //     os << "VTK"
                //        << "/" << VtkBaseName << ".0000" << imgm[RealNo] << ".vtk" << ends;
                // else if (imgm[RealNo] < 100)
                //     os << "VTK"
                //        << "/" << VtkBaseName << ".000" << imgm[RealNo] << ".vtk" << ends;
                // else if (imgm[RealNo] < 1000)
                //     os << "VTK"
                //        << "/" << VtkBaseName << ".00" << imgm[RealNo] << ".vtk" << ends;
                // else if (imgm[RealNo] < 10000)
                //     os << "VTK"
                //        << "/" << VtkBaseName << ".0" << imgm[RealNo] << ".vtk" << ends;
                // else
                //     os << "VTK"
                //        << "/" << VtkBaseName << "." << imgm[RealNo] << ".vtk" << ends;
                // Output->WriteVtk(os.str().c_str());

                // }
                // if (imgm[RealNo] < 10)
                //     fileoutMC = "MonteCarlo/" + filename + "_t0000" + std::to_string(imgm[RealNo]) + ".txt";
                // else if (imgm[RealNo] < 100)
                //     fileoutMC = "MonteCarlo/" + filename + "_t000" + std::to_string(imgm[RealNo]) + ".txt";
                // else if (imgm[RealNo] < 1000)
                //     fileoutMC = "MonteCarlo/" + filename + "_t00" + std::to_string(imgm[RealNo]) + ".txt";
                // else if (imgm[RealNo] < 10000)
                //     fileoutMC = "MonteCarlo/" + filename + "_t0" + std::to_string(imgm[RealNo]) + ".txt";
                // else
                //     fileoutMC = "MonteCarlo/" + filename + "_t" + std::to_string(imgm[RealNo]) + ".txt";
                // imgm[RealNo]++;
                // printToTxt(fileoutMC, sol, N_DOF, 1, 'C');
                fileoutMC = generateFileName("MonteCarlo/" + filename, imgTxt[RealNo], N_Realisations);
                printToTxt(fileoutMC, sol, N_DOF, 1, 'C');
                imgTxt[RealNo]++;

            } // for(l=0;l< N_SubSteps;l++)

        } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

        //======================================================================
        // produce final outout
        //======================================================================
        printVTKOutput(VtkBaseName, imgVtk + RealNo, Output);

        // os.seekp(std::ios::beg);
        // if (imgm[RealNo] < 10)
        //     os << "VTK"
        //        << "/" << VtkBaseName << ".0000" << imgm[RealNo] << ".vtk" << ends;
        // else if (imgm[RealNo] < 100)
        //     os << "VTK"
        //        << "/" << VtkBaseName << ".000" << imgm[RealNo] << ".vtk" << ends;
        // else if (imgm[RealNo] < 1000)
        //     os << "VTK"
        //        << "/" << VtkBaseName << ".00" << imgm[RealNo] << ".vtk" << ends;
        // else if (imgm[RealNo] < 10000)
        //     os << "VTK"
        //        << "/" << VtkBaseName << ".0" << imgm[RealNo] << ".vtk" << ends;
        // else
        //     os << "VTK"
        //        << "/" << VtkBaseName << "." << imgm[RealNo] << ".vtk" << ends;
        // Output->WriteVtk(os.str().c_str());
        // if (imgm[RealNo] < 10)
        //     fileoutMC = "MonteCarlo/" + filename + "_t0000" + std::to_string(imgm[RealNo]) + ".txt";
        // else if (imgm[RealNo] < 100)
        //     fileoutMC = "MonteCarlo/" + filename + "_t000" + std::to_string(imgm[RealNo]) + ".txt";
        // else if (imgm[RealNo] < 1000)
        //     fileoutMC = "MonteCarlo/" + filename + "_t00" + std::to_string(imgm[RealNo]) + ".txt";
        // else if (imgm[RealNo] < 10000)
        //     fileoutMC = "MonteCarlo/" + filename + "_t0" + std::to_string(imgm[RealNo]) + ".txt";
        // else
        //     fileoutMC = "MonteCarlo/" + filename + "_t" + std::to_string(imgm[RealNo]) + ".txt";
        // imgm[RealNo]++;
        // printToTxt(fileoutMC, sol, N_DOF, 1, 'C');
        fileoutMC = generateFileName("MonteCarlo/" + filename, imgTxt[RealNo], N_Realisations);
        printToTxt(fileoutMC, sol, N_DOF, 1, 'C');
        imgTxt[RealNo]++;
        // cout << " Solution Norm After: " << Ddot(N_DOF,sol,sol) <<endl;

        // set Current Time as Zero
        TDatabase::TimeDB->CURRENTTIME = 0;
        // delete SystemMatrix;
    }
    delete[] CompositeVectorMC;
    delete[] imgVtk;
    delete[] imgTxt;
    return;
}