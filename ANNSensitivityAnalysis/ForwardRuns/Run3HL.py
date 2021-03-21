import numpy as np
import os
from dataFunctions import *

def runSimulations(thisRunDir):
# Activation function array for output layer
# Note: 0: sigmoid, 1:identity, 2:leakyReLU, 3:tanH
    OPLTYPE_ARRAY = [0,1,3,4];

# Dimension array for the hidden layer
    HL_DIM_ARRAY = [5,7,10];

# Activation function array for the hidden layer
# Note: 0: sigmoid, 1:identity, 2:ReLU, 3:leakyReLU, 4:tanH
    HL_TYPE_ARRAY = [ 0,3,4];

# Epochs array
    EPOCHS_ARRAY = [10,100];

# Number of hidden layers fixed at 3
    NHL = 3;


# Initialization a dictionary of all the parameters
    data = dict();
    data["NHL"] = NHL;
    data["OPLTYPE"] = 1;
    data["HL_0_DIM"] = 5;
    data["HL_0_TYPE"] = 3;
    data["HL_1_DIM"] = 5;
    data["HL_1_TYPE"] = 3;
    data["HL_2_DIM"] = 5;
    data["HL_2_TYPE"] = 3;
    data["IPDATADIM"] = 1;
    data["IPLDIM"] = 3;
    data["OPLDIM"] = 1;
    data["DATASET_NAME"] = "trainingData.csv";
    data["TRAINING_DATA_PERCENTAGE"] = 73;
    data["VALIDATION_DATA_PERCENTAGE"] = 7;
    data["OPTIMIZER_CODE"] = 3;
    data["OPTIMIZER_STEP_SIZE"] = 0.001;
    data["SGD_BATCH_SIZE"] = 32;
    data["MAX_ITERATIONS"] = 1000000;
    data["TOLERANCE"] = 1e-8;
    data["DROPOUT_RATIO"] = 0.2;
    data["EPOCHS"] = 100;
    data["SAVE_DATA_FILE"] = "testResults.csv";



# Start the loops

    listDir = os.listdir(thisRunDir);
    maxDir = 0;
    if (listDir == ['metadata.dat', 'inputSpace.dat']):
        maxDir = 0;
    else:
        listDir.remove('metadata.dat');
        listDir.remove('inputSpace.dat');
        listDir.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
        maxDir = int(listDir[-1]);
        maxDir = maxDir + 1;

    count = maxDir;
    for OPLTYPE in OPLTYPE_ARRAY:
        for HL_0_DIM in HL_DIM_ARRAY:
            for HL_0_TYPE in HL_TYPE_ARRAY:
                for HL_1_DIM in HL_DIM_ARRAY:
                    for HL_1_TYPE in HL_TYPE_ARRAY:
                        for HL_2_DIM in HL_DIM_ARRAY:
                            for HL_2_TYPE in HL_TYPE_ARRAY:
                                for EPOCHS in EPOCHS_ARRAY:
                                    # Create a directory to save the configuration and the result
                                    
                                    os.mkdir(thisRunDir+'/' +str(count));

                                    # Update the dictionary for this run
                                    data["OPLTYPE"] = OPLTYPE;
                                    data["HL_0_DIM"] = HL_0_DIM;
                                    data["HL_0_TYPE"] = HL_0_TYPE;
                                    data["HL_1_DIM"] = HL_1_DIM;
                                    data["HL_1_TYPE"] = HL_1_TYPE;
                                    data["HL_2_DIM"] = HL_2_DIM;
                                    data["HL_2_TYPE"] = HL_2_TYPE;
                                    data["EPOCHS"] = EPOCHS;
                                    data["SAVE_DATA_FILE"] = "testResults.csv";

                                    # Create the code script with this data 
                                    createScript(data);

                                    # Send a copy of this script to the respective results directory
                                    os.system("cp runScript.dat "+thisRunDir+'/'+str(count)+'/');

                                    # Write the input parameter space specifications in inputSapceFile
                                    inputSpaceFile = open(thisRunDir+'/inputSpace.dat', 'a');
                                    inputSpaceFile.write(str(count)+'   '+  str(NHL)+'   '+ str(OPLTYPE)+ '   ' + str(HL_0_DIM) + '   ' + str(HL_0_TYPE) +  '   ' + str(HL_1_DIM) + '   ' + str(HL_1_TYPE) + '   ' + str(HL_2_DIM) + '   '+ str(HL_2_TYPE) + '   '+  str(EPOCHS) + '   \n');

                                    inputSpaceFile.close();
                                    # Run the simulation
                                    os.system('./parmoon_2D_SEQUENTIAL.exe runScript.dat');

                                    # Copy the results file to the respective folder
                                    os.system("mv ' testResults.csv' "+thisRunDir+'/'+str(count)+'/testResults.csv');

                                    print(count);

                                    # Update the file
                                    count += 1;


