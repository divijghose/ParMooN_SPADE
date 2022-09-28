/**
 * @brief Routine to print data to a text file for post-processing
 *
 * @param filename Name of output file
 * @param printArray Matrix to be printed
 * @param height Number of rows in printArray
 * @param width Number of columns in printArray
 * @param RowOrColMaj 'R' if data is stored in Row Major form, 'C' if column major
 */
void printToTxt(std::string filename, double *printArray, int height, int width, char RowOrColMaj)
{
    std::ofstream printFile;
    printFile.open(filename);
    if (RowOrColMaj == 'R')
    {
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                printFile << printArray[i * width + j];
                if (j != width - 1)
                    printFile << ",";
            }
            printFile << endl;
        }
    }
    else if (RowOrColMaj == 'C')
    {
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                printFile << printArray[j * height + i];
                if (j != width - 1)
                    printFile << ",";
            }
            printFile << endl;
        }
    }
    printFile.close();
    cout << "File printed succesfully: " << filename << endl;
}

/**
 * @brief Routine to read the Monte Carlo realization matrix from a text file
 *
 * @param RealznVect Pointer to the realization matrix
 * @param N_R Number of Realizations
 * @param N_DOF Number of degrees of freedom
 */
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
        cout << "Could not open the file\n";

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

/**
 * @brief Routine to write the realization matrix to a text file
 *
 * @param RealznVect Pointer to realization matrix
 * @param N_R Number of realizations
 * @param N_DOF Number of degrees of freedom
 */
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

/**
 * @brief Routine to generate a formatted file name
 *
 * @param baseName Base name of file
 * @param m Time step for which file is is being written
 * @param N_R Number of realizations for which the file is being written
 * @return std::string
 */
std::string generateFileName(std::string baseName, int m, int N_R)
{
    std::string fileName;
    if (m < 10)
        fileName = baseName + std::to_string(N_R) + "_t0000" + std::to_string(m) + ".txt";
    else if (m < 100)
        fileName = baseName + std::to_string(N_R) + "_t000" + std::to_string(m) + ".txt";
    else if (m < 1000)
        fileName = baseName + std::to_string(N_R) + "_t00" + std::to_string(m) + ".txt";
    else if (m < 10000)
        fileName = baseName + std::to_string(N_R) + "_t0" + std::to_string(m) + ".txt";
    else
        fileName = baseName + std::to_string(N_R) + "_t" + std::to_string(m) + ".txt";

    return fileName;
}

/**
 * @brief Routine to read data from a text file
 * 
 * @param fileName Name of the text file
 * @param Vector Pointer to store matrix
 * @param height Number of rows in the matrix
 * @param width Number of columns in the matrix
 * @param RowOrColMaj 'R' if matrix is stored in row major, 'C' if column major
 */
void readFromText(std::string fileName, double *Vector, int height, int width, char RowOrColMaj)
{
    cout << "Read In" << endl;
    std::vector<std::vector<std::string>> content;
    std::vector<std::string> row;
    std::string line, word;

    std::ifstream file(fileName);
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
        cout << "File " << fileName << " opened succesfully" << endl;
    }
    else
        cout << "Could not open the file " << fileName << "\n";

    cout << "" << endl;
    if (RowOrColMaj == 'R')
    {
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                Vector[i * width + j] = std::stod(content[i][j]);
            }
        }
    }
    else if (RowOrColMaj == 'C')
    {
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                Vector[j * height + i] = std::stod(content[i][j]);
            }
        }
    }

    cout << "File " << fileName << " read successfully" << endl;
    return;
}

