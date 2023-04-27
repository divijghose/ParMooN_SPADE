
void GenerateRealizations(TFESpace2D *Velocity_FeSpace, TFESpace2D *Pressure_FeSpace, double *RealizationVector)
{
    // ///////////////////////////////////////////////////////////////////////////////////////////////
    // ////////// -------- REALISATION DATA GENERATION ----------------------------------------- //////
    // ///////////////////////////////////////////////////////////////////////////////////////////////

    int N_Realisations = TDatabase::ParamDB->REALIZATIONS;
    int N_DOF = Velocity_FeSpace->GetN_DegreesOfFreedom();
    int Nx, Ny; // Number of grid points in x and y directions

    Nx = sqrt(N_DOF);
    Ny = sqrt(N_DOF);
    cout << " Nx : " << Nx << " Ny : " << Ny << endl;

    double dx = 1.0 / (Nx - 1);
    double dy = 1.0 / (Ny - 1);
    double width, height;

    if (TDatabase::ParamDB->toggleRealznSource == 0)
    {
        double LengthScale = TDatabase::ParamDB->LENGTHSCALE;
        double EigenPercent = TDatabase::ParamDB->EIGENPERCENT;

        double *x_coord_true, *y_coord_true, *x_coord_calc, *y_coord_calc;

        x_coord_true = new double[N_DOF]();
        y_coord_true = new double[N_DOF]();
        x_coord_calc = new double[N_DOF]();
        y_coord_calc = new double[N_DOF]();

        int *mappingArray = new int[N_DOF]();

        int N = (TDatabase::ParamDB->ANSATZ_ORDER * pow(2, TDatabase::ParamDB->UNIFORM_STEPS)) + 1;
        // int N = pow(2, TDatabase::ParamDB->UNIFORM_STEPS) + 1;

        cout << " N : " << N << endl;

        for (int i = 0; i < N_DOF; i++)
        {
            int local_i = i / N;
            int local_j = i % N;

            x_coord_calc[i] = double(1.0 / (N - 1)) * local_i;
            y_coord_calc[i] = double(1.0 / (N - 1)) * local_j;
        }

        Velocity_FeSpace->GetDOFPosition(x_coord_true, y_coord_true);

        for (int i = 0; i < N_DOF; i++) // Calculated Values
        {
            bool foundFlag = false;

            for (int j = 0; j < N_DOF; j++) // True Values
            {
                if (abs(x_coord_calc[i] - x_coord_true[j]) < 1e-10 && abs(y_coord_calc[i] - y_coord_true[j]) < 1e-10)
                {
                    mappingArray[i] = j;
                    foundFlag = true;
                }
            }

            if (!foundFlag)
            {
                cerr << " DOF NOT FOUND FOR " << i << " position : " << setw(8) << x_coord_true[i] << setw(8) << y_coord_true[i] << endl;
                cout << "Exiting...Please check order of elements and uniform steps" << endl;
                exit(0);
            }
        }
        OutPut("True values match calculated values for all co-ordinates" << endl);

        // if x_coord_true[i] and y_coord_true[i] are the true coordinates of the ith DOF, then mappingArray[i] is the index of the ith DOF in the calculated coordinates

        double *x_coord_int, *y_coord_int;
        x_coord_int = new double[(Nx - 2) * (Ny - 2)]();
        y_coord_int = new double[(Nx - 2) * (Ny - 2)]();

        // fill x_coord_int with all values of x_coord_calc except the boundary values
        // fill y_coord_int with all values of y_coord_calc except the boundary values

        int int_count = 0;
        int *mappingArrayInternal = new int[(Nx - 2) * (Ny - 2)]();
        for (int i = 0; i < N_DOF; i++)
        {
            if (x_coord_calc[i] != 0 && x_coord_calc[i] != 1 && y_coord_calc[i] != 0 && y_coord_calc[i] != 1)
            {
                x_coord_int[int_count] = x_coord_calc[i];
                y_coord_int[int_count] = y_coord_calc[i];
                int_count++;
            }
        }

        cout << "int_count : " << int_count << endl;
        cout << "Nx_int :" << sqrt(int_count) << " Ny_int : " << sqrt(int_count) << endl;

        for (int i = 0; i < int_count; i++) // Calculated Values
        {
            bool foundFlag = false;

            for (int j = 0; j < N_DOF; j++) // True Values
            {
                if (abs(x_coord_int[i] - x_coord_true[j]) < 1e-10 && abs(y_coord_int[i] - y_coord_true[j]) < 1e-10)
                {
                    mappingArrayInternal[i] = j;
                    foundFlag = true;
                }
            }

            if (!foundFlag)
            {
                cerr << " Internal DOF NOT FOUND FOR " << i << " position : " << setw(8) << x_coord_true[i] << setw(8) << y_coord_true[i] << endl;
                cout << "Exiting...Please check order of elements and uniform steps" << endl;
                exit(0);
            }
        }
        OutPut("True values match calculated values for all co-ordinates of Internal DOFs" << endl);

        // ///////////////////////////////////////////////////////////////////////////////////////////////
        // ////////// -------- MOLLIFIER FUNCTION (To have consistent boundary conditions)   ------- //////
        // ///////////////////////////////////////////////////////////////////////////////////////////////

        int kmax = 5; // Number of averaging steps - Chamge to higher value for smoother realizations.

        double *wgt = new double[(Nx - 2) * (Ny - 2)]();                    // Mollifier function
        double *wgt_int = new double[(Nx - 2 - 2) * (Ny - 2 - 2)]();    // Mollifier function on interior window // Row major
        double *wgt_right = new double[(Nx - 2 - 2) * (Ny - 2 - 2)]();  // Mollifier function on right shift window // Row major
        double *wgt_left = new double[(Nx - 2 - 2) * (Ny - 2 - 2)]();   // Mollifier function on left shift window // Row major
        double *wgt_top = new double[(Nx - 2 - 2) * (Ny - 2 - 2)]();    // Mollifier function on top shift window // Row major
        double *wgt_bottom = new double[(Nx - 2 - 2) * (Ny - 2 - 2)](); // Mollifier function on bottom shift window    // Row major
        double *wgt_rt = new double[(Nx - 2 - 2) * (Ny - 2 - 2)]();     // Molliifier function on right top shift window // Row major
        double *wgt_rb = new double[(Nx - 2 - 2) * (Ny - 2 - 2)]();     // Molliifier function on right bottom shift window     // Row major
        double *wgt_lt = new double[(Nx - 2 - 2) * (Ny - 2 - 2)]();     // Molliifier function on left top shift window // Row major
        double *wgt_lb = new double[(Nx - 2 - 2) * (Ny - 2 - 2)]();     // Molliifier function on left bottom shift window // Row major

        for (int a = 0; a < int_count; a++)
        {
            if (1 - x_coord_int[a] <= 1e-6 || 1 - y_coord_int[a] <= 1e-6 || x_coord_int[a] - 0 <= 1e-6 || y_coord_int[a] - 0 <= 1e-6)
                wgt[a] = 0; // 0 on the boundary
            else
                wgt[a] = 1; // 1 on the interior
        }

        for (int k = 0; k < kmax; k++) // Averaging pass
        {
            for (int b = 0; b < Ny - 2 - 2; b++)
            {
                for (int d = 0; d < Nx - 2 - 2; d++)
                {
                    wgt_int[b * (Nx - 2 - 2) + d] = wgt[((b + 1) * Nx - 2) + (d + 1)];
                    wgt_right[b * (Nx - 2 - 2) + d] = wgt[((b + 1) * Nx - 2) + (d + 2)];
                    wgt_left[b * (Nx - 2 - 2) + d] = wgt[((b + 1) * Nx - 2) + (d)];
                    wgt_top[b * (Nx - 2 - 2) + d] = wgt[((b)*Nx - 2) + (d + 1)];
                    wgt_bottom[b * (Nx - 2 - 2) + d] = wgt[((b + 2) * Nx - 2) + (d + 1)];
                    wgt_rt[b * (Nx - 2 - 2) + d] = wgt[((b)*Nx - 2) + (d + 2)];
                    wgt_rb[b * (Nx - 2 - 2) + d] = wgt[((b + 2) * Nx - 2) + (d + 2)];
                    wgt_lt[b * (Nx - 2 - 2) + d] = wgt[((b)*Nx - 2) + (d)];
                    wgt_lb[b * (Nx - 2 - 2) + d] = wgt[((b + 2) * Nx - 2) + (d)];
                }
            }
            for (int b = 0; b < Ny - 2 - 2; b++)
            {
                for (int d = 0; d < Nx - 2 - 2; d++)
                {
                    wgt[((b + 1) * Nx - 2) + (d + 1)] *= (wgt_right[b * (Nx - 2 - 2) + d] + wgt_left[b * (Nx - 2 - 2) + d] + wgt_top[b * (Nx - 2 - 2) + d] + wgt_bottom[b * (Nx - 2 - 2) + d] + wgt_rt[b * (Nx - 2 - 2) + d] + wgt_rb[b * (Nx - 2 - 2) + d] + wgt_lt[b * (Nx - 2 - 2) + d] + wgt_lb[b * (Nx - 2 - 2) + d]) / 8;
                }
            }
        }
        printToTxt("wgt.txt", wgt, Ny - 2, Nx - 2, 'R');
        
        ////////////////////////////////////////////////////////////////////////////////////
        ////////////////// Correlation Matrix Calculation //////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////
        double *CorrMatx = new double[N_DOF * N_DOF](); // Correlation Matrix
        double *CovMatx = new double[N_DOF * N_DOF]();  // Covariance Matrix

        double r = 0.0;
        double norm = 0;

        for (int i = 0; i < N_DOF; i++)
        {
            double x_1_fix = x_coord_int[i];
            double y_1_fix = y_coord_int[i];

            for (int j = 0; j < N_DOF; j++)
            {
                double x_2_move = x_coord_int[j];
                double y_2_move = y_coord_int[j];

                r = sqrt(pow((x_1_fix - x_2_move), 2) + pow((y_1_fix - y_2_move), 2)); // Compute pair-wise distances

                // Correlation Calculation

                CorrMatx[i * N_DOF + j] = exp((-1.0 * r) / (LengthScale)) * (1 + (r / LengthScale) + pow(r / (LengthScale), 2) / 3.0) * 1e-6;

                // Covariance Calculation

                CovMatx[i * N_DOF + j] = wgt[i] * wgt[j] * CorrMatx[i * N_DOF + j];
            }
        }
        std::string fileOutCorrelation = "Init/Correlation_" + std::to_string(N_DOF) + "_NR_" + std::to_string(N_Realisations) + ".txt";
        printToTxt(fileOutCorrelation, CorrMatx, N_DOF, N_DOF, 'R');

        std::string fileOutCovariance = "Init/Covariance_" + std::to_string(N_DOF) + "_NR_" + std::to_string(N_Realisations) + ".txt";
        printToTxt(fileOutCovariance, CovMatx, N_DOF, N_DOF, 'R');

        ////////////////////////////////////////////////////// SVD ////////////////////////////////////////////
        // Declare SVD parameters
        MKL_INT m1 = N_DOF, n = N_DOF, lda = N_DOF, ldu = N_DOF, ldvt = N_DOF, info;
        double superb[std::min(N_DOF, N_DOF) - 1];
        double w[N];

        double *S = new double[N_DOF]();
        double *U = new double[N_DOF * N_DOF]();
        double *Vt = new double[N_DOF * N_DOF]();

        // info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', N_DOF, C, N_DOF, S);

        info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', m1, n, CovMatx, lda,
                              S, U, ldu, Vt, ldvt, superb);

        if (info > 0)
        {
            printf("The algorithm computing SVD of Covariance Matrix failed to converge.\n");
            exit(1);
        }

        std::string fileOutCorrS = "Init/S_of_Cov" + std::to_string(N_DOF) + "_NR_" + std::to_string(N_Realisations) + ".txt";
        printToTxt(fileOutCorrS, S, N_DOF, 1, 'C');

        std::string fileOutCorrU = "Init/U_of_Cov" + std::to_string(N_DOF) + "_NR_" + std::to_string(N_Realisations) + ".txt";
        printToTxt(fileOutCorrU, U, N_DOF, N_DOF, 'R');

        std::string fileOutCorrVt = "Init/Vt_of_Cov" + std::to_string(N_DOF) + "_NR_" + std::to_string(N_Realisations) + ".txt";
        printToTxt(fileOutCorrVt, Vt, N_DOF, N_DOF, 'R');

        int energyVal = 0;
        int temp = 0;

        double sumSingularVal = 0;
        for (int i = 0; i < N_DOF; i++)
            sumSingularVal += S[i];

        double val = 0;
        for (energyVal = 0; energyVal < N_DOF; energyVal++)
        {
            val += S[energyVal];
            temp++;
            if (val / sumSingularVal > 0.99)
                break;
        }

        int modDim = temp + 1;
        cout << " MODES : " << modDim << endl;

        double *Ut = new double[N_DOF * modDim]();
        double *Z = new double[N_Realisations * modDim]();

        double *RealizationVectorTemp = new double[N_DOF * N_Realisations]();

        // -------------- Generate Random Number Based on Normal Distribution -------------------------//
        int k = 0;
        int skip = N_DOF - modDim;
        int count = 0;
        for (int i = 0; i < N_DOF * N_DOF; i++)
        {
            if (count < modDim)
            {
                Ut[k] = U[i];

                count++;
                k++;
            }
            else
            {
                i += skip;
                count = 0;
                i--;
            }
        }

        double varNormFactor = S[0];
        // for (int k = 0; k < modDim; k++)
        //     S[k] = (S[k] / varNormFactor);

        for (int k = 0; k < modDim; k++)
        {
            std::random_device rd{};
            std::mt19937 gen{rd()};
            std::normal_distribution<> d{0, sqrt(S[k])};

            double *norm1 = new double[N_Realisations]();

            for (int n = 0; n < N_Realisations; ++n)
            {
                Z[k * N_Realisations + n] = d(gen);
            }
        }

        cout << " N_Realisations : " << N_Realisations << endl;
        cout << " MULT START " << endl;
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_DOF, N_Realisations, modDim, 1.0, Ut, modDim, Z, N_Realisations, 0.0, RealizationVector, N_Realisations);
        cout << " MULT DONE " << endl;

        for (int i = 0; i < N_DOF; i++)
        {
            for (int j = 0; j < N_Realisations; j++)
            {
                RealizationVectorTemp[mappingArray[i] * N_Realisations + j] = RealizationVector[j + N_Realisations * i];
            }
        }

        memcpy(RealizationVector, RealizationVectorTemp, N_DOF * N_Realisations * SizeOfDouble);

        cout << N_Realisations << " REALISATIONS COMPUTED " << endl;

        delete[] Ut;
        delete[] Z;
        delete[] RealizationVectorTemp;
        delete[] S;
        delete[] U;
        delete[] Vt;
        delete[] x_coord_true;
        delete[] y_coord_true;
        delete[] x_coord_calc;
        delete[] y_coord_calc;
        delete[] CorrMatx;
        delete[] CovMatx;

        if (TDatabase::ParamDB->writeRealznToText == 1)
            writeRealizationToText(RealizationVector, N_Realisations, N_DOF);
    }
    else if (TDatabase::ParamDB->toggleRealznSource == 1)
        readRealizationFromText(RealizationVector, N_Realisations, N_DOF);
    else
    {
        cout << "Please select correct value of Realization_Source" << endl
             << TDatabase::ParamDB->toggleRealznSource << " is not an acceptable value" << endl;
        exit(0);
    }

    /////////////////////////////////////// -------- END OF REALISATION DATA SETS ------------ ////////////////////////////////////////////////////////////////
}
