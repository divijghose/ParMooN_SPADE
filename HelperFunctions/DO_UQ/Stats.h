
void calcMeanRealization(const double *RealznVect, double *MeanVect, const int N_R, const int N_DOF)
{
    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = 0; j < N_R; j++)
        {
            MeanVect[i] += RealznVect[i * N_R + j] / N_R;
        }
    }

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

    return;
}


double calc_MeanFieldEnergy(TFESpace2D *Fespace, TFEFunction2D *FeScalar_Cmean)
{
    // double val = 0;
    double mfe = 0;
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

    int lenMean = FeScalar_Cmean->GetLength();
    double *C_Mean_Array = FeScalar_Cmean->GetValues();

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

        // Save Values of C at all quadrature points for I component
        double C_Mean[N_Points2];

        for (int quadPt = 0; quadPt < N_Points2; quadPt++)
            C_Mean[quadPt] = 0;

        for (int quadPt = 0; quadPt < N_Points2; quadPt++)
        {
            for (int j = 0; j < N_BaseFunct; j++)
            {
                int globDOF = DOF[j];
                C_Mean[quadPt] += origvaluesD00[quadPt][j] * C_Mean_Array[globDOF];
            }
        }

        // INner Quadrature Loop
        for (int quadPt = 0; quadPt < N_Points2; quadPt++)
        {
            double Mult = Weights2[quadPt] * AbsDetjk[quadPt];
            mfe += C_Mean[quadPt] * C_Mean[quadPt] * Mult; // This is Final "f"

        } // inner quadrature loop

        // --
    } // cell loop

    return mfe;
} // calc_MFE function end
void calc_ModeOrtho(TFESpace2D *Fespace, TFEVectFunct2D *FeVector_Cmode, int N_S, double *IPMatxFE)
{
    // double val = 0;
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

    double *C_Mode_Array_i = new double[lenMode * N_S]();
    double *C_Mode_Array_j = new double[lenMode * N_S]();

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
                    for (int x = 0; x < N_BaseFunct; x++)
                    {
                        int globDOF = DOF[x];
                        C_j[quadPt] += origvaluesD00[quadPt][x] * C_Mode_Array_j[globDOF];
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
    }     // cell loop
    delete[] C_Mode_Array_i;
    delete[] C_Mode_Array_j;
    return;
} // calc_MFE function end

void calc_princVariance(double *princVariance, int N_S)
{

    double *localCov = new double[N_S * N_S]();
    memcpy(localCov, TDatabase::ParamDB->COVARIANCE_MATRIX_DO, N_S * N_S * SizeOfDouble);
    double *wi, *VL, *VR;
    int info, lwork;
    double wkopt;
    double *work;
    wi = new double[N_S]();
    VL = new double[N_S * N_S]();
    VL = new double[N_S * N_S]();
    VR = new double[N_S * N_S]();
    int N = N_S;
    int LDA = N_S;
    int LDVL = N_S;
    int LDVR = N_S;
    double *wr = new double[N_S]();
    lwork = -1;

    info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'N', N_S, localCov, LDA, princVariance, wi,
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
}


void reconstructMCfromDO(double *recon, double *meanDO, double *coeffDO, double *modeDO, int N_R, int N_DOF, int N_S)
{
    double *PertVect = new double[N_DOF * N_R]();
    double *modeRowMaj = new double[N_DOF * N_S]();
    for (int i = 0; i < N_DOF; i++)
    {
        for (int j = 0; j < N_S; j++)
        {
            modeRowMaj[i * N_S + j] = modeDO[j * N_DOF + i];
        }
    }
    double *coeffRowMaj = new double[N_R * N_S]();
    for (int i = 0; i < N_R; i++)
    {
        for (int j = 0; j < N_S; j++)
        {
            coeffRowMaj[i * N_S + j] = coeffDO[j * N_R + i];
        }
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N_DOF, N_R, N_S, 1.0, modeRowMaj, N_S, coeffRowMaj, N_S, 0.0, PertVect, N_R);
    for (int i = 0; i < N_DOF; ++i)
    {
        for (int j = 0; j < N_R; ++j)
        {
            recon[i * N_R + j] = PertVect[i * N_R + j] + meanDO[i];
        }
    }
    delete[] PertVect;
    return;
}


