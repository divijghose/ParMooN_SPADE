
void GenerateRealizations(TFESpace2D *Scalar_FeSpace, double *RealizationVector)
{
    // ///////////////////////////////////////////////////////////////////////////////////////////////
    // ////////// -------- REALISATION DATA GENERATION ----------------------------------------- //////
    // ///////////////////////////////////////////////////////////////////////////////////////////////

    int N_Realisations = TDatabase::ParamDB->REALIZATIONS;
    int N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
    int i;
    if (TDatabase::ParamDB->toggleRealznSource == 0)
    {
        double LengthScale = TDatabase::ParamDB->LENGTHSCALE;
        double EigenPercent = TDatabase::ParamDB->EIGENPERCENT;

        double *org_x_coord = new double[N_DOF]();
        double *org_y_coord = new double[N_DOF]();
        double *x_coord = new double[N_DOF]();
        double *y_coord = new double[N_DOF]();
        int *mappingArray = new int[N_DOF]();

        i = 0;
  
        int N = (TDatabase::ParamDB->ANSATZ_ORDER*pow(2, TDatabase::ParamDB->UNIFORM_STEPS)) + 1;
        for (int i = 0; i < N_DOF; i++)
        {
            int local_i = i / N;
            int local_j = i % N;

            x_coord[i] = double(1.0 / (N - 1)) * local_i;
            y_coord[i] = double(1.0 / (N - 1)) * local_j;
        }
        Scalar_FeSpace->GetDOFPosition(org_x_coord, org_y_coord);

        for (int i = 0; i < N_DOF; i++) // Generated Values
        {
            // get the generated Value
            double xx = x_coord[i];
            double yy = y_coord[i];
            bool foundFlag = false;

            for (int j = 0; j < N_DOF; j++) // Actual parmooN Co-ordinates
            {
                if (abs(xx - org_x_coord[j]) < 1e-10 && abs(yy - org_y_coord[j]) < 1e-10)
                {
                    mappingArray[i] = j;
                    foundFlag = true;
                }
            }

            if (!foundFlag)
                cerr << " DOF NOT FOUND FOR " << i << " position : " << setw(8) << org_x_coord[i] << setw(8) << org_y_coord[i] << endl;
        }
        // int N_DOF =  N * N;
        double *x = new double[N_DOF]();
        double *y = new double[N_DOF]();

        for (int i = 0; i < N_DOF; i++)
        {
            int local_i = i / N;
            int local_j = i % N;

            x[i] = double(1.0 / (N - 1)) * local_j;
            y[i] = double(1.0 / (N - 1)) * local_i;
        }
        int kmax = 4;
        int Nx, Ny;
        Nx = sqrt(N_DOF);
        Ny = sqrt(N_DOF);
        double *wgt = new double[Nx * Ny]();
        double *wgt_int = new double[(Nx - 2) * (Ny - 2)]();
        double *wgt_right = new double[(Nx - 2) * (Ny - 2)]();
        double *wgt_left = new double[(Nx - 2) * (Ny - 2)]();
        double *wgt_top = new double[(Nx - 2) * (Ny - 2)]();
        double *wgt_bottom = new double[(Nx - 2) * (Ny - 2)]();
        double *wgt_rt = new double[(Nx - 2) * (Ny - 2)]();
        double *wgt_rb = new double[(Nx - 2) * (Ny - 2)]();
        double *wgt_lt = new double[(Nx - 2) * (Ny - 2)]();
        double *wgt_lb = new double[(Nx - 2) * (Ny - 2)]();
        if (TDatabase::ParamDB->stddev_switch == 5)
        {

            for (int a = 0; a < N_DOF; a++)
            {
                if (1 - x[a] <= 1e-6 || 1 - y[a] <= 1e-6 || x[a] - 0 <= 1e-6 || y[a] - 0 <= 1e-6)
                    wgt[a] = 0;
                else
                    wgt[a] = 1;
            }

            for (int k = 0; k < kmax; k++)
            {
                for (int b = 0; b < Ny - 2; b++)
                {
                    for (int d = 0; d < Nx - 2; d++)
                    {
                        wgt_int[b * (Nx - 2) + d] = wgt[((b + 1) * Nx) + (d + 1)];
                        wgt_right[b * (Nx - 2) + d] = wgt[((b + 1) * Nx) + (d + 2)];
                        wgt_left[b * (Nx - 2) + d] = wgt[((b + 1) * Nx) + (d)];
                        wgt_top[b * (Nx - 2) + d] = wgt[((b)*Nx) + (d + 1)];
                        wgt_bottom[b * (Nx - 2) + d] = wgt[((b + 2) * Nx) + (d + 1)];
                        wgt_rt[b * (Nx - 2) + d] = wgt[((b)*Nx) + (d + 2)];
                        wgt_rb[b * (Nx - 2) + d] = wgt[((b + 2) * Nx) + (d + 2)];
                        wgt_lt[b * (Nx - 2) + d] = wgt[((b)*Nx) + (d)];
                        wgt_lb[b * (Nx - 2) + d] = wgt[((b + 2) * Nx) + (d)];
                    }
                }
                for (int b = 0; b < Ny - 2; b++)
                {
                    for (int d = 0; d < Nx - 2; d++)
                    {
                        wgt[((b + 1) * Nx) + (d + 1)] *= (wgt_right[b * (Nx - 2) + d] + wgt_left[b * (Nx - 2) + d] + wgt_top[b * (Nx - 2) + d] + wgt_bottom[b * (Nx - 2) + d] + wgt_rt[b * (Nx - 2) + d] + wgt_rb[b * (Nx - 2) + d] + wgt_lt[b * (Nx - 2) + d] + wgt_lb[b * (Nx - 2) + d]) / 8;
                    }
                }
            }
            printToTxt("wgt.txt", wgt, Ny, Nx, 'R');
        }
        else if (TDatabase::ParamDB->stddev_switch == 6)
        {
            double disp = TDatabase::ParamDB->stddev_disp;
            double halfwidth = 0.2;
            double xcenter, ycenter;
            xcenter = 0.0;
            ycenter = 0.0;
            double xdist, ydist, refxdist, refydist;
            refxdist = 2.0;
            refydist = 2.0;
            for (int m = 0; m < N_DOF; m++)
            {
                xdist = fabs(x[m] - disp);
                if (xdist <= refxdist)
                {
                    xcenter = x[m];
                    refxdist = xdist;
                }
            }
            for (int m = 0; m < N_DOF; m++)
            {
                ydist = fabs(y[m] - disp);
                if (ydist <= refydist)
                {
                    ycenter = y[m];
                    refydist = ydist;
                }
            }

            int Nx, Ny;
            Nx = sqrt(N_DOF);
            Ny = sqrt(N_DOF);

            for (int k = 0; k < N_DOF; k++)
            {
                if (fabs(x[k] - xcenter) <= halfwidth && fabs(y[k] - ycenter) <= halfwidth)
                    wgt[k] = 1;
            }

            kmax = 6;
            for (int k = 0; k < kmax; k++)
            {
                for (int b = 0; b < Ny - 2; b++)
                {
                    for (int d = 0; d < Nx - 2; d++)
                    {
                        wgt_int[b * (Nx - 2) + d] = wgt[((b + 1) * Nx) + (d + 1)];
                        wgt_right[b * (Nx - 2) + d] = wgt[((b + 1) * Nx) + (d + 2)];
                        wgt_left[b * (Nx - 2) + d] = wgt[((b + 1) * Nx) + (d)];
                        wgt_top[b * (Nx - 2) + d] = wgt[((b)*Nx) + (d + 1)];
                        wgt_bottom[b * (Nx - 2) + d] = wgt[((b + 2) * Nx) + (d + 1)];
                        wgt_rt[b * (Nx - 2) + d] = wgt[((b)*Nx) + (d + 2)];
                        wgt_rb[b * (Nx - 2) + d] = wgt[((b + 2) * Nx) + (d + 2)];
                        wgt_lt[b * (Nx - 2) + d] = wgt[((b)*Nx) + (d)];
                        wgt_lb[b * (Nx - 2) + d] = wgt[((b + 2) * Nx) + (d)];
                    }
                }
                for (int b = 0; b < Ny - 2; b++)
                {
                    for (int d = 0; d < Nx - 2; d++)
                    {
                        wgt[((b + 1) * Nx) + (d + 1)] *= (wgt_right[b * (Nx - 2) + d] + wgt_left[b * (Nx - 2) + d] + wgt_top[b * (Nx - 2) + d] + wgt_bottom[b * (Nx - 2) + d] + wgt_rt[b * (Nx - 2) + d] + wgt_rb[b * (Nx - 2) + d] + wgt_lt[b * (Nx - 2) + d] + wgt_lb[b * (Nx - 2) + d]) / 8;
                    }
                }
            }
            printToTxt("wgt.txt", wgt, Ny, Nx, 'R');
        }
        else if (TDatabase::ParamDB->stddev_switch == 7)
        {
            double disp = TDatabase::ParamDB->stddev_disp;
            double halfwidth = 0.2;
            double xcenter, ycenter;
            xcenter = 0.0;
            ycenter = 0.0;
            double xdist, ydist, refxdist, refydist;
            refxdist = 2.0;
            refydist = 2.0;
            for (int m = 0; m < N_DOF; m++)
            {
                xdist = fabs(x[m] - disp);
                if (xdist <= refxdist)
                {
                    xcenter = x[m];
                    refxdist = xdist;
                }
            }
            for (int m = 0; m < N_DOF; m++)
            {
                ydist = fabs(y[m] - disp);
                if (ydist <= refydist)
                {
                    ycenter = y[m];
                    refydist = ydist;
                }
            }

            int Nx, Ny;
            Nx = sqrt(N_DOF);
            Ny = sqrt(N_DOF);

            for (int k = 0; k < N_DOF; k++)
            {
                if (sqrt(pow(x[k] - xcenter, 2) + pow(y[k] - ycenter, 2)) <= halfwidth)
                    wgt[k] = 1;
            }

            kmax = 7;
            for (int k = 0; k < kmax; k++)
            {
                for (int b = 0; b < Ny - 2; b++)
                {
                    for (int d = 0; d < Nx - 2; d++)
                    {
                        wgt_int[b * (Nx - 2) + d] = wgt[((b + 1) * Nx) + (d + 1)];
                        wgt_right[b * (Nx - 2) + d] = wgt[((b + 1) * Nx) + (d + 2)];
                        wgt_left[b * (Nx - 2) + d] = wgt[((b + 1) * Nx) + (d)];
                        wgt_top[b * (Nx - 2) + d] = wgt[((b)*Nx) + (d + 1)];
                        wgt_bottom[b * (Nx - 2) + d] = wgt[((b + 2) * Nx) + (d + 1)];
                        wgt_rt[b * (Nx - 2) + d] = wgt[((b)*Nx) + (d + 2)];
                        wgt_rb[b * (Nx - 2) + d] = wgt[((b + 2) * Nx) + (d + 2)];
                        wgt_lt[b * (Nx - 2) + d] = wgt[((b)*Nx) + (d)];
                        wgt_lb[b * (Nx - 2) + d] = wgt[((b + 2) * Nx) + (d)];
                    }
                }
                for (int b = 0; b < Ny - 2; b++)
                {
                    for (int d = 0; d < Nx - 2; d++)
                    {
                        wgt[((b + 1) * Nx) + (d + 1)] *= (wgt_right[b * (Nx - 2) + d] + wgt_left[b * (Nx - 2) + d] + wgt_top[b * (Nx - 2) + d] + wgt_bottom[b * (Nx - 2) + d] + wgt_rt[b * (Nx - 2) + d] + wgt_rb[b * (Nx - 2) + d] + wgt_lt[b * (Nx - 2) + d] + wgt_lb[b * (Nx - 2) + d]) / 8;
                    }
                }
            }
            printToTxt("wgt.txt", wgt, Ny, Nx, 'R');
        }

        double *C = new double[N_DOF * N_DOF]();  // MATRIX
        double *C1 = new double[N_DOF * N_DOF](); // MATRIX  - Corelation Matrix

        double r = 0.0;
        double norm = 0;

        for (int i = 0; i < N_DOF; i++)
        {
            double actual_x = x[i];
            double actual_y = y[i];

            for (int j = 0; j < N_DOF; j++)
            {
                double local_x = x[j];
                double local_y = y[j];

                r = sqrt(pow((actual_x - local_x), 2) + pow((actual_y - local_y), 2));

                // CO -Relation
                C[i * N_DOF + j] = exp((-1.0 * r * 10.0) / (LengthScale)) * (1.0 + (10.0 * r / LengthScale) + (pow(10.0 * r / (LengthScale), 2) / 3.0));
                C1[i * N_DOF + j] = exp((-1.0 * r) / (LengthScale)) * (1 + (r / LengthScale) + pow(r / (LengthScale), 2) / 3.0);

                if (TDatabase::ParamDB->stddev_switch == 0)
                {
                    double sig_r1 = exp(-1.0 / (1.0 - pow((2 * actual_x - 1), 4))) * exp(-1.0 / (1 - pow((2 * actual_y - 1), 4)));
                    double sig_r2 = exp(-1.0 / (1.0 - pow((2 * local_x - 1), 4))) * exp(-1.0 / (1 - pow((2 * local_y - 1), 4)));
                    // Co Variance
                    C[i * N_DOF + j] *= sig_r1 * sig_r2 * 5.0;
                }

                else if (TDatabase::ParamDB->stddev_switch == 1)
                {

                    double E = TDatabase::ParamDB->stddev_denom;
                    double disp = TDatabase::ParamDB->stddev_disp;
                    double power = TDatabase::ParamDB->stddev_power;
                    double sig_r1 = (exp(-1.0 * pow((2 * actual_x - 1 - disp), power) / (E)) / (2 * Pi * sqrt(E))) * (exp(-1.0 * pow((2 * actual_y - 1 - disp), power) / (E)) / (2 * Pi * sqrt(E)));
                    double sig_r2 = (exp(-1.0 * pow((2 * local_x - 1 - disp), power) / (E)) / (2 * Pi * sqrt(E))) * (exp(-1.0 * pow((2 * local_y - 1 - disp), power) / (E)) / (2 * Pi * sqrt(E)));
                    // Co Variance
                    C[i * N_DOF + j] *= 1.0 * sig_r1 * sig_r2;
                }

                else if (TDatabase::ParamDB->stddev_switch == 2)
                {
                    double amplitude = TDatabase::ParamDB->stddev_power;
                    double sig_r1 = (amplitude)*sin(-1.0 * Pi * (2 * actual_x - 2)) * sin(-1.0 * Pi * (2 * actual_y - 2));
                    double sig_r2 = (amplitude)*sin(-1.0 * Pi * (2 * local_x - 2)) * sin(-1.0 * Pi * (2 * local_y - 2));
                    C[i * N_DOF + j] *= sig_r1 * sig_r2;
                }

                else if (TDatabase::ParamDB->stddev_switch == 3)
                {

                    double sigma = 1 / (LengthScale * (2.0 * sqrt(Pi)));
                    double E = TDatabase::ParamDB->stddev_denom;
                    double disp = TDatabase::ParamDB->stddev_disp;
                    double power = TDatabase::ParamDB->stddev_power;
                    double sig_r1 = sigma * (exp(-1.0 * pow((2 * actual_x - 1 - disp), power) / (E * pow(LengthScale, 2)))) * (exp(-1.0 * pow((2 * actual_y - 1 - disp), power) / (E * pow(LengthScale, 2))));
                    double sig_r2 = sigma * (exp(-1.0 * pow((2 * local_x - 1 - disp), power) / (E * pow(LengthScale, 2)))) * (exp(-1.0 * pow((2 * local_y - 1 - disp), power) / (E * pow(LengthScale, 2))));
                    // Co Variance
                    C[i * N_DOF + j] *= 1 * sig_r1 * sig_r2;
                }
                else if (TDatabase::ParamDB->stddev_switch == 4)
                {

                    double height = 1;
                    double disp = TDatabase::ParamDB->stddev_disp;
                    double spread = TDatabase::ParamDB->stddev_denom;
                    double power = TDatabase::ParamDB->stddev_power;
                    double sig_r1 = height * (exp(-1.0 * (((pow(2.0 * actual_x - 1.0 - disp, power)) / (pow(spread, 2))))) * (exp(-1.0 * ((pow(2.0 * actual_y - 1.0 - disp, power)) / (pow(spread, 2))))));
                    double sig_r2 = height * (exp(-1.0 * (((pow(2.0 * local_x - 1.0 - disp, power)) / (pow(spread, 2))))) * (exp(-1.0 * ((pow(2 * local_y - 1.0 - disp, power)) / (pow(spread, 2))))));
                    // Co Variance
                    C[i * N_DOF + j] *= sig_r1 * sig_r2 * 5e-3;
                }
                else if (TDatabase::ParamDB->stddev_switch == 5)
                {

                    double sig_r1, sig_r2;
                    sig_r1 = wgt[i];
                    sig_r2 = wgt[j];
                    // Co Variance
                    // C[i * N_DOF + j] *= sig_r1 * sig_r2;
                    C[i * N_DOF + j] *= sig_r1 * sig_r2;
                }
                else if (TDatabase::ParamDB->stddev_switch == 6)
                {
                    double sig_r1, sig_r2;
                    sig_r1 = wgt[i];
                    sig_r2 = wgt[j];
                    // Co Variance
                    // C[i * N_DOF + j] *= sig_r1 * sig_r2;
                    C[i * N_DOF + j] *= sig_r1 * sig_r2;
                }
                else if (TDatabase::ParamDB->stddev_switch == 7)
                {
                    double sig_r1, sig_r2;
                    sig_r1 = wgt[i];
                    sig_r2 = wgt[j];
                    // Co Variance
                    // C[i * N_DOF + j] *= sig_r1 * sig_r2;
                    C[i * N_DOF + j] *= sig_r1 * sig_r2*700;
                }
                else
                {
                    cout << "Error - No standard deviation function is defined for stddev_switch: " << TDatabase::ParamDB->stddev_switch << endl;
                    exit(0);
                }

                // norm += C[j * N + i] * C[j * N + i];
            }
        }

        std::string fileOutCorrelation = "Init/Correlation_" + std::to_string(N_DOF) + "_NR_" + std::to_string(N_Realisations) + ".txt";
        printToTxt(fileOutCorrelation, C, N_DOF, N_DOF, 'R');

        ////////////////////////////////////////////////////// SVD ////////////////////////////////////////////
        // Declare SVD parameters
        MKL_INT m1 = N_DOF, n = N_DOF, lda = N_DOF, ldu = N_DOF, ldvt = N_DOF, info;
        double superb[std::min(N_DOF, N_DOF) - 1];
        double w[N];

        double *S = new double[N_DOF]();
        double *U = new double[N_DOF * N_DOF]();
        double *Vt = new double[N_DOF * N_DOF]();

        // info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', N_DOF, C, N_DOF, S);

        info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', m1, n, C, lda,
                              S, U, ldu, Vt, ldvt, superb);

        if (info > 0)
        {
            printf("The algorithm computing SVD of Covariance Matrix failed to converge.\n");
            exit(1);
        }

        std::string fileOutCorrS = "Init/Corr_S_" + std::to_string(N_DOF) + "_NR_" + std::to_string(N_Realisations) + ".txt";
        printToTxt(fileOutCorrS, S, N_DOF, 1, 'C');

        std::string fileOutCorrU = "Init/Corr_U_" + std::to_string(N_DOF) + "_NR_" + std::to_string(N_Realisations) + ".txt";
        printToTxt(fileOutCorrU, U, N_DOF, N_DOF, 'R');

        std::string fileOutCorrVt = "Init/Corr_Vt_" + std::to_string(N_DOF) + "_NR_" + std::to_string(N_Realisations) + ".txt";
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

        cout << " MODES : " << temp + 1 << endl;

        int modDim = temp + 1;

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
        delete[] org_x_coord;
        delete[] org_y_coord;
        delete[] x_coord;
        delete[] y_coord;
        delete[] C;
        delete[] C1;

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