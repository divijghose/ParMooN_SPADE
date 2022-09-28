
/**
 * @brief Routine to claculate the inner product of a matrix
 * 
 * @param IPMatx Pointer to the inner product 
 * @param Vector Pointer to the matrix
 * @param height Number of rows of matrix
 * @param width Number of columns of matrix
 * @param RoworColMaj 'R' if matrix is row major, 'C' if column major
 */
void calcIPMatx(double *IPMatx, double *Vector, int height, int width, char RoworColMaj)
{
    if (RoworColMaj == 'R')
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, width, width, height, 1, Vector, width, Vector, width, 0.0, IPMatx, width);
    else if (RoworColMaj == 'C')
    {
        double *TempRowMaj = new double[height * width]();
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                TempRowMaj[i * width + j] = Vector[j * height + i];
            } // end for j
        }     // end for i
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, width, width, height, 1, TempRowMaj, width, TempRowMaj, width, 0.0, IPMatx, width);
    }
}


void normalizeStochasticModes(TFESpace2D *Fespace, TFEVectFunct2D *FeVector_Cmode, int N_S, double *C_Stoch_Norm)
{

    double *IPMatxFE = new double[N_S * N_S]();
    int N_Cells = Fespace->GetN_Cells();
    TCollection *coll = Fespace->GetCollection();

    // Get the Global DOF arrays INdex from the FE Space.
    int *GlobalNumbers = Fespace->GetGlobalNumbers();
    int *BeginIndex = Fespace->GetBeginIndex();

    // --- Quadrature Formula Arrays  ------------------//
    int N_Points2;
    double *Weights2, *t1, *t2; // Weights - Quadrature Weights array ; t1  -- Quadrature point ( xi ) in ref coordinate ; t2 -  Quadrature Point ( eta ) in Ref Coordinate
    bool Needs2ndDer[1];
    Needs2ndDer[0] = TRUE;
    double AbsDetjk[MaxN_PointsForNodal2D];
    double X[MaxN_PointsForNodal2D];
    double Y[MaxN_PointsForNodal2D];

    // FE Values Arrays
    double **origvaluesD00; // Shape function values at quadrature Points
    double **origvaluesD10; // Shape Function Derivatives ( x ) at Quadrature Points
    double **origvaluesD01; // Shape Function Derivatives ( y ) at Quadrature Points
    double **origvaluesD20; // Shape Function 2nd Derivatives ( x ) at Quadrature Points
    double **origvaluesD02; // Shape Function 2nd Derivatives ( y ) at Quadrature Points

    int lenMode = FeVector_Cmode->GetLength();
    double *C_Mode_Array = FeVector_Cmode->GetValues();

    double *C_Mode_Array_i = new double[lenMode]();
    double *C_Mode_Array_j = new double[lenMode]();

    double *stochFactors = new double[lenMode * N_S]();

    for (int cellId = 0; cellId < N_Cells; cellId++)
    { // cell loop
        TBaseCell *currentCell = coll->GetCell(cellId);
        // Get the "ID" of Finite Element for the given 2D Element ( Conforming/NonConforming-Order Finite Element : eg : it could be Conforming-2nd order Finite Element )
        FE2D elementId = Fespace->GetFE2D(cellId, currentCell);
        // Get the Class object for that 2d FEM Element , which has all the details like Shape functions , Nodal point locations for that location, Reference Transformation ( Affine , Bilinear )
        TFE2D *element = TFEDatabase2D::GetFE2D(elementId);
        TFEDesc2D *fedesc = element->GetFEDesc2D();
        // Class for basis functions in 2D ( Number of basis functions ), basis function values and Derivatives
        TBaseFunct2D *bf = element->GetBaseFunct2D();
        // Get the Reference Elemet
        BF2DRefElements RefElement = TFEDatabase2D::GetRefElementFromFE2D(elementId);
        // Get the reference Transformation -- Affine Mapping / Bilnea Mapping of Triangle or Quadrilateral
        RefTrans2D referenceTransformation = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(elementId);
        // Get the number of basis functions in the Current Cell ( Number of Local DOF)
        int N_BaseFunct = element->GetN_DOF();
        // Type of Basis Function in 2D
        BaseFunct2D BaseFunct_ID = element->GetBaseFunct2D_ID();

        // get cell measure
        double hK = currentCell->GetDiameter();

        switch (referenceTransformation)
        {
        case QuadBilinear:
        {
            int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
            QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
            TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
            QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

            // Set the values on the Reference Cell
            TRefTrans2D *F_K = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
            TFEDatabase2D::SetCellForRefTrans(currentCell, QuadBilinear); // Set the Cell for Current reference Transformation

            // Get Original Coordinates from reference Coordinates and the Determinant of jacobian
            TFEDatabase2D::GetOrigFromRef(QuadBilinear, N_Points2, t1, t2, X, Y, AbsDetjk); // Get the Original Co-orinates for the cell from xi values

            // Get all the original Values from the Referece cell values.
            TFEDatabase2D::GetOrigValues(QuadBilinear, 1, &BaseFunct_ID, N_Points2, t1, t2, QF2, Needs2ndDer);

            // The below are 2D arrays in the form
            // Values[QuadraturePointLocation][ShapeFunction]  i.e, the Value of Shapefunction at all quadrature points for each shape functions
            origvaluesD00 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D00); // Shape Function Values at Quadrature Points
            origvaluesD10 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D10); // Shape Function Derivative Values at Quadrature Points
            origvaluesD01 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D01); // Shape Function Derivative Values at Quadrature Point
            origvaluesD20 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D20); // Shape Function 2nd Derivative Values at Quadrature Point
            origvaluesD02 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D02); // Shape Function 2nd Derivative Values at Quadrature Point
            break;
        }

        case QuadAffin:
        {
            int l = bf->GetPolynomialDegree();                                     // Get the Polynomial Degreee  of the basis functions
            QuadFormula2D QF2 = TFEDatabase2D::GetQFQuadFromDegree(3 * l);         // Get te ID of Quadrature Formula
            TQuadFormula2D *QuadratureRule = TFEDatabase2D::GetQuadFormula2D(QF2); // Get the Quadrature Rule Objetc based on Quadrature ID
            QuadratureRule->GetFormulaData(N_Points2, Weights2, t1, t2);           // get the Quadrature points , Weights

            // Set the values on the Reference Cell
            TRefTrans2D *F_K = TFEDatabase2D::GetRefTrans2D(QuadAffin);
            TFEDatabase2D::SetCellForRefTrans(currentCell, QuadAffin); // Set the Cell for Current reference Transformation

            // Get Original Coordinates from reference Coordinates and the Determinant of jacobian
            TFEDatabase2D::GetOrigFromRef(QuadAffin, N_Points2, t1, t2, X, Y, AbsDetjk); // Get the Original Co-orinates for the cell from xi values

            // Get all the original Values from the Referece cell values.
            TFEDatabase2D::GetOrigValues(QuadAffin, 1, &BaseFunct_ID, N_Points2, t1, t2, QF2, Needs2ndDer);

            // The below are 2D arrays in the form
            // Values[QuadraturePointLocation][ShapeFunction]  i.e, the Value of Shapefunction at all quadrature points for each shape functions
            origvaluesD00 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D00); // Shape Function Values at Quadrature Points
            origvaluesD10 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D10); // Shape Function Derivative Values at Quadrature Points
            origvaluesD01 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D01); // Shape Function Derivative Values at Quadrature Point
            origvaluesD20 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D20); // Shape Function 2nd Derivative Values at Quadrature Points
            origvaluesD02 = TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D02); // Shape Function 2nd Derivative Values at Quadrature Point

            break;
        }

        default:
        {
            cout << " [ERROR] - Error in File : CoeffEqn_DO.C " << endl;
            cout << " Unknown Reftype " << endl;
            cout << " REF TYPE : " << referenceTransformation << endl;
            exit(0);
            break;
        }
        }

        int *DOF = GlobalNumbers + BeginIndex[cellId];

        for (int i = 0; i < N_S; i++)
        { // i-loop start
            memcpy(C_Mode_Array_i, C_Mode_Array + i * lenMode, lenMode * SizeOfDouble);
            double C_i[N_Points2];

            for (int quadPt = 0; quadPt < N_Points2; quadPt++)
                C_i[quadPt] = 0.;
            // Obtain all values for C_a

            for (int quadPt = 0; quadPt < N_Points2; quadPt++)
            { // quadpt loop start
                for (int nb = 0; nb < N_BaseFunct; nb++)
                {
                    int globDOF = DOF[nb];
                    C_i[quadPt] += origvaluesD00[quadPt][nb] * C_Mode_Array_i[globDOF];
                    stochFactors[globDOF + i * lenMode] += 1;
                }
            } // quadpt loop end

            for (int j = 0; j < N_S; j++)
            { // j-loop start

                memcpy(C_Mode_Array_j, C_Mode_Array + j * lenMode, lenMode * SizeOfDouble);
                double C_j[N_Points2];
                for (int quadPt = 0; quadPt < N_Points2; quadPt++)
                    C_j[quadPt] = 0.;
                // Obtain all values for C_a
                for (int quadPt = 0; quadPt < N_Points2; quadPt++)
                { // quadpt loop start
                    for (int nb = 0; nb < N_BaseFunct; nb++)
                    {
                        int globDOF = DOF[nb];
                        C_j[quadPt] += origvaluesD00[quadPt][nb] * C_Mode_Array_j[globDOF];
                    }
                } // quadpt loop end

                // INner Quadrature Loop
                for (int quadPt = 0; quadPt < N_Points2; quadPt++)
                {
                    double Mult = Weights2[quadPt] * AbsDetjk[quadPt];
                    IPMatxFE[i + j * N_S] += C_i[quadPt] * C_j[quadPt] * Mult;

                } // Inner Quadrature Loop

            } // j-loop end

        } // i-loop end

        // --
    } // cell loop

    for (int j = 0; j < N_S; j++)
    {
        for (int i = 0; i < lenMode; i++)
        {
            C_Stoch_Norm[i + j * lenMode] = C_Mode_Array[i + j * lenMode] / sqrt(IPMatxFE[j + j * N_S]);
        }
    }
    double *IPNew = new double[N_S * N_S]();
    calcIPMatx(IPNew, C_Stoch_Norm, lenMode, N_S, 'C');
    for (int j = 0; j < N_S; j++)
    {
        for (int i = 0; i < lenMode; i++)
        {
            C_Stoch_Norm[i + j * lenMode] = C_Stoch_Norm[i + j * lenMode] / sqrt(IPNew[j + j * N_S]);
        }
    }
    calcIPMatx(IPNew, C_Stoch_Norm, lenMode, N_S, 'C');
    // for (int j = 0; j < N_S; j++)
    // {
    //     for (int i = 0; i < lenMode; i++)
    //     {
    //         C_Stoch_Norm[i + j * lenMode] = C_Mode_Array[i + j * lenMode];
    //     }
    // }

    printToTxt("stoch.txt", IPNew, N_S, N_S, 'C');
    delete[] C_Mode_Array_i;
    delete[] C_Mode_Array_j;

    delete[] IPMatxFE;

    return;
}


void qr(double *const _Q, double *const _R, double *const _A, const size_t _m, const size_t _n)
{
    // Maximal rank is used by Lapacke
    const size_t rank = std::min(_m, _n);

    // Tmp Array for Lapacke
    const std::unique_ptr<double[]> tau(new double[rank]);

    // Calculate QR factorisations
    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, (int)_m, (int)_n, _A, (int)_n, tau.get());

    // Copy the upper triangular Matrix R (rank x _n) into position
    for (size_t row = 0; row < rank; ++row)
    {
        memset(_R + row * _n, 0, row * sizeof(double));                                // Set starting zeros
        memcpy(_R + row * _n + row, _A + row * _n + row, (_n - row) * sizeof(double)); // Copy upper triangular part from Lapack result.
    }

    // Create orthogonal matrix Q (in tmpA)
    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, (int)_m, (int)rank, (int)rank, _A, (int)_n, tau.get());

    // Copy Q (_m x rank) into position
    if (_m == _n)
    {
        memcpy(_Q, _A, sizeof(double) * (_m * _n));
    }
    else
    {
        for (size_t row = 0; row < _m; ++row)
        {
            memcpy(_Q + row * rank, _A + row * _n, sizeof(double) * (rank));
        }
    }
    return;
}



void reorthonormalizeA(double *Mode, int N_DOF, int N_S)
{

    double *Temp = new double[N_DOF * N_S]();
    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = 0; j < N_S; j++)
        {
            Temp[i * N_S + j] = Mode[j * N_DOF + i];
        }
    }
    double *Q = new double[N_DOF * N_S]();
    double *R = new double[N_DOF * N_S]();
    qr(Q, R, Temp, N_DOF, N_S);
    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = 0; j < N_S; j++)
        {
            Mode[j * N_DOF + i] = Q[i * N_S + j];
        }
    }
    delete[] Temp;
    delete[] Q;
    delete[] R;

    return;
}

void reorthonormalizeB(double *Mode, double *Coeff, int N_DOF, int N_S, int N_R)
{
    cout << "Enter new reorthonormalization routine" << endl;
    int height = N_R;
    int width = N_S;
    double *Cov = new double[N_S * N_S]();
    double *phi = new double[width * height](); // Col to Row Major
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            phi[width * i + j] = Coeff[i + height * j];
        }
    }

    const double divVal = (1.0 / (height - 1));
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, width, width, height, divVal, phi, width, phi, width, 0.0, Cov, width);
    cout << "Covariance for new reorthonormalization routine calculated" << endl;

    double *wi, *VDCL, *VDCR;
    int info, lwork;
    double wkopt;
    double *work;
    wi = new double[N_S]();
    VDCL = new double[N_S * N_S]();
    VDCR = new double[N_S * N_S]();
    int N = N_S;
    int LDA = N_S;
    int LDVL = N_S;
    int LDVR = N_S;
    double *wr = new double[N_S]();
    lwork = -1;
    double *DDC = new double[N_S]();

    info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'V', 'N', N_S, Cov, LDA, DDC, wi,
                         VDCL, LDVL, VDCR, LDVR);

    if (info == 0)
        cout << "The routine computing eignevalues of coefficient matrix was successful" << endl;
    else if (info < 0)
    {
        cout << "The routine computing eignevalues of coefficient matrix was unsuccessful" << endl
             << -1 * info << "th argument invalid" << endl;
        exit(0);
    }
    else if (info > 0)
    {
        cout << "the QR algorithm failed to compute all the eigenvalues" << endl;

        exit(0);
    }

    double *ModeDC = new double[N_DOF * N_S]();
    double *CoeffDC = new double[N_R * N_S]();
    double *ModeRow = new double[N_DOF * N_S]();
    double *CoeffRow = new double[N_R * N_S]();
    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = 0; j < N_S; j++)
        {
            ModeRow[i * N_S + j] = Mode[j * N_DOF + i];
        }
    }
    for (int i = 0; i < N_R; i++)
    {
        for (int j = 0; j < N_S; j++)
        {
            CoeffRow[i * N_S + j] = Coeff[j * N_R + i];
        }
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_DOF, N_S, N_S, 1.0, ModeRow, N_S, VDCL, N_S, 0.0, ModeDC, N_S);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_R, N_S, N_S, 1.0, CoeffRow, N_S, VDCL, N_S, 0.0, CoeffDC, N_S);

    double *M = new double[N_S * N_S]();
    calcIPMatx(M, ModeDC, N_DOF, N_S, 'R');

    double *VML = new double[N_S * N_S]();
    double *VMR = new double[N_S * N_S]();
    double *DM = new double[N_S]();

    info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'V', 'N', N_S, M, LDA, DM, wi,
                         VML, LDVL, VMR, LDVR);

    if (info == 0)
        cout << "The routine computing eignevalues of coefficient matrix was successful" << endl;
    else if (info < 0)
    {
        cout << "The routine computing eignevalues of coefficient matrix was unsuccessful" << endl
             << -1 * info << "th argument invalid" << endl;
        exit(0);
    }
    else if (info > 0)
    {
        cout << "the QR algorithm failed to compute all the eigenvalues" << endl;

        exit(0);
    }

    double *CoeffOC = new double[N_R * N_S]();
    double *ModeOC = new double[N_DOF * N_S]();

    double *D = new double[N_S * N_S]();
    double *Dinv = new double[N_S * N_S]();

    for (int i = 0; i < N_S; i++)
    {
        D[i * N_S + i] = DM[i];
        Dinv[i * N_S + i] = DM[i];
    }

    matInv(Dinv, N_S);
    for (int i = 0; i < N_S; i++)
    {
        D[i * N_S + i] = sqrt(D[i * N_S + i]);
        Dinv[i * N_S + i] = sqrt(Dinv[i * N_S + i]);
    }

    double *Temp = new double[N_S * N_S]();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_S, N_S, N_S, 1.0, VML, N_S, D, N_S, 0.0, Temp, N_S);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_R, N_S, N_S, 1.0, CoeffDC, N_S, Temp, N_S, 0.0, CoeffOC, N_S);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_S, N_S, N_S, 1.0, VML, N_S, Dinv, N_S, 0.0, Temp, N_S);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_DOF, N_S, N_S, 1.0, ModeDC, N_S, Temp, N_S, 0.0, ModeOC, N_S);

    double *CovOC = new double[N_S * N_S]();
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, width, width, height, divVal, CoeffOC, width, CoeffOC, width, 0.0, CovOC, width);

    double *VOL = new double[N_S * N_S]();
    double *VOR = new double[N_S * N_S]();
    double *DDO = new double[N_S]();

    info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'V', 'N', N_S, CovOC, LDA, DDO, wi,
                         VOL, LDVL, VOR, LDVR);

    if (info == 0)
        cout << "The routine computing eignevalues of coefficient matrix was successful" << endl;
    else if (info < 0)
    {
        cout << "The routine computing eignevalues of coefficient matrix was unsuccessful" << endl
             << -1 * info << "th argument invalid" << endl;
        exit(0);
    }
    else if (info > 0)
    {
        cout << "the QR algorithm failed to compute all the eigenvalues" << endl;

        exit(0);
    }
    double tracedc = 0.0;
    double tracedo = 0.0;
    for (int i = 0; i < N_S; i++)
    {
        tracedc += DDC[i];
        tracedo += DDO[i];
    }
    double Mult = sqrt(tracedc / tracedo);
    double *CoeffO = new double[N_R * N_S]();
    double *ModeO = new double[N_DOF * N_S]();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_R, N_S, N_S, Mult, CoeffOC, N_S, VOL, N_S, 0.0, CoeffO, N_S);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_DOF, N_S, N_S, 1.0, ModeOC, N_S, VOL, N_S, 0.0, ModeO, N_S);

    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = 0; j < N_S; j++)
        {
            Mode[j * N_DOF + i] = ModeO[i * N_S + j];
        }
    }
    for (int i = 0; i < N_S; i++)
    {
        for (int j = 0; j < N_S; j++)
        {
            Coeff[j * N_S + i] = CoeffO[i * N_S + j] * Mult;
        }
    }

    delete[] M;
    delete[] Cov;
    delete[] phi;
    delete[] wi;
    delete[] VDCL;
    delete[] VDCR;
    delete[] wr;
    delete[] DDC;
    delete[] ModeRow;
    delete[] CoeffDC;
    delete[] ModeDC;
    delete[] D;
    delete[] Dinv;
    delete[] Temp;
    delete[] CoeffOC;
    delete[] CoeffO;
    delete[] ModeOC;
    delete[] ModeO;
    delete[] VOL;
    delete[] VOR;
    delete[] DDO;
    delete[] CovOC;

    return;
}

void reorthonormalizeC(double *Mode, int N_DOF, int N_S)
{
    double *Temp = new double[N_S * N_S]();
    double *P = new double[N_S * N_S]();

    calcIPMatx(Temp, Mode, N_DOF, N_S, 'C');
    for (int i = 0; i < N_S; i++)
    {
        for (int j = 0; j < N_S; j++)
        {
            P[i * N_S + j] = Temp[j * N_S + i];
        }
    }

    double *wi, *VL, *VR;
    int info, lwork;
    double wkopt;
    double *work;
    wi = new double[N_S]();
    VL = new double[N_S * N_S]();
    VR = new double[N_S * N_S]();
    int N = N_S;
    int LDA = N_S;
    int LDVL = N_S;
    int LDVR = N_S;
    double *wr = new double[N_S]();
    lwork = -1;
    double *D = new double[N_S]();

    info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'V', 'V', N_S, P, LDA, D, wi,
                         VL, LDVL, VR, LDVR);

    if (info == 0)
        cout << "The routine computing eignevalues of coefficient matrix was successful" << endl;
    else if (info < 0)
    {
        cout << "The routine computing eignevalues of coefficient matrix was unsuccessful" << endl
             << -1 * info << "th argument invalid" << endl;
        exit(0);
    }
    else if (info > 0)
    {
        cout << "the QR algorithm failed to compute all the eigenvalues" << endl;

        exit(0);
    }

    double *RootD = new double[N_S * N_S]();
    for (int i = 0; i < N_S; i++)
        RootD[i + i * N_S] = sqrt(D[i]);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_S, N_S, N_S, 1.0, VL, N_S, RootD, N_S, 0.0, Temp, N_S);

    double *RootP = new double[N_S * N_S]();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N_S, N_S, N_S, 1.0, Temp, N_S, VL, N_S, 0.0, RootP, N_S);

    matInv(RootP, N_S);

    double *ModeOld = new double[N_DOF * N_S]();
    double *ModeNew = new double[N_DOF * N_S]();

    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = 0; j < N_S; j++)
        {
            ModeOld[i * N_S + j] = Mode[j * N_DOF + i];
        }
    }

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_DOF, N_S, N_S, 1.0, ModeOld, N_S, RootP, N_S, 0.0, ModeNew, N_S);

    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = 0; j < N_S; j++)
        {
            Mode[j * N_DOF + i] = ModeNew[i * N_S + j];
        }
    }
    delete[] wi;
    delete[] wr;
    delete[] D;
    delete[] P;
    delete[] Temp;
    delete[] VR;
    delete[] VL;
    delete[] ModeNew;
    delete[] ModeOld;
    delete[] RootD;
    delete[] RootP;

    return;
}

