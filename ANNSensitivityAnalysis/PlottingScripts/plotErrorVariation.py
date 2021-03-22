import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dataFunctions as DF

from plotError import plotError

def getError(projectName, runNumber):
    curDir = os.getcwd();
    #_______________________________________________________
    # Set paths and variables
    #_______________________________________________________
    # Output path for sensitivity analysis (./output)
    outputPath = os.getcwd()+'/output';
    # Location for storing test case (name is specified in inputData.py)
    projectOutputDir = outputPath+'/'+projectName; 
    runDir = projectOutputDir+'/'+str(runNumber); 

    # Change into the run directory inside output/caseName
    os.chdir(runDir);

    #_______________________________________________________
    # Prepare Sample data
    #_______________________________________________________

    #DF.createOutputSpace(runDir);

    # Read metadata
    # Load input space for sensitivity analysis
    inputData = np.loadtxt("inputSpace.dat");
    # Total number of samples
    numberOfSamples = inputData.shape[0];

    # Load the output space for sensitivity analysis
    outputData = np.loadtxt("outputSpace.dat");

    #_______________________________________________________
    # Process the data 
    #_______________________________________________________
    p5 = np.percentile(outputData,5, axis=0);
    p5={'L1Error':p5[0], 'L2Error':p5[1], 'MinError':p5[2], 'MaxError':p5[3], 'MSError':p5[4]};

    p95 = np.percentile(outputData,95, axis=0);
    p95={'L1Error':p95[0], 'L2Error':p95[1], 'MinError':p95[2], 'MaxError':p95[3], 'MSError':p95[4]};

    pMean = np.mean(outputData, axis=0);
    pMean={'L1Error':pMean[0], 'L2Error':pMean[1], 'MinError':pMean[2], 'MaxError':pMean[3], 'MSError':pMean[4]};

    pStd = np.std(outputData, axis=0);
    pStd={'L1Error':pStd[0], 'L2Error':pStd[1], 'MinError':pStd[2], 'MaxError':pStd[3], 'MSError':pStd[4]};

    os.chdir(curDir);
    return (pMean, p5, p95, pStd);


if __name__ == "__main__":

    curDir = os.getcwd();

    TotalRuns = 7;

    # 4 metrics: mean value, 5 percentile, 95 percentile, standard deviation
    L1Error = np.zeros(shape=(TotalRuns, 4));
    L2Error = np.zeros(shape=(TotalRuns, 4));
    MSError = np.zeros(shape=(TotalRuns, 4));
    MinError = np.zeros(shape=(TotalRuns, 4));
    MaxError = np.zeros(shape=(TotalRuns, 4));

    projectName = "Expt1";

    projectDir = os.getcwd()+"/output/"+projectName+"/";

    for runNumber in range(TotalRuns):
        print (runNumber);
        (pMean, p5, p95, pStd) = getError(projectName , runNumber);

        L1Error[runNumber, 0] = pMean['L1Error'];
        L2Error[runNumber, 0] = pMean['L2Error'];
        MinError[runNumber, 0] = pMean['MinError'];
        MaxError[runNumber, 0] = pMean['MaxError'];
        MSError[runNumber, 0] = pMean['MSError'];

        L1Error[runNumber,  1] = p5['L1Error'];
        L2Error[runNumber,  1] = p5['L2Error'];
        MinError[runNumber, 1] = p5['MinError'];
        MaxError[runNumber, 1] = p5['MaxError'];
        MSError[runNumber,  1] = p5['MSError'];

        L1Error[runNumber,  2] = p95['L1Error'];
        L2Error[runNumber,  2] = p95['L2Error'];
        MinError[runNumber, 2] = p95['MinError'];
        MaxError[runNumber, 2] = p95['MaxError'];
        MSError[runNumber,  2] = p95['MSError'];

        L1Error[runNumber,  3] = pStd['L1Error'];
        L2Error[runNumber,  3] = pStd['L2Error'];
        MinError[runNumber, 3] = pStd['MinError'];
        MaxError[runNumber, 3] = pStd['MaxError'];
        MSError[runNumber,  3] = pStd['MSError'];

        pass;

    #_______________________________________________________
    # plot the data: 
    #_______________________________________________________

    x = np.arange(TotalRuns);

    fig, axs = plt.subplots(2,2, figsize=(6.4,4.8), dpi=300, constrained_layout=True);

    # Title:
    plt.suptitle("Variation of error vs. size of training dataset");

    #_______________________________________________________
    # L1 error plot
    #_______________________________________________________
    location = (0,0);
    axs[location].plot(x, L1Error[:,0], '-', color ='black', linewidth='1.5', label = 'Mean error');
    axs[location].plot(x, L1Error[:,1], '--', color ='green', linewidth='2', label = '5-Percentile error');
    axs[location].plot(x, L1Error[:,2], '-.', color ='red', linewidth='2', label = '95-Percentile error');
    axs[location].plot(x, L1Error[:,3], '.', color ='blue', linewidth='3', label = 'Standard deviation');

    axs[location].set_xticks([0,1,2,3,4,5,6]);
    axs[location].set_xticklabels(['12','25','50','100','200','400','800']);

    axs[location].set_ylabel(r"$L_1$ Error");
    axs[location].legend(loc='upper left', fontsize=5);

    #_______________________________________________________
    # MS error plot
    #_______________________________________________________
    location = (1,0);
    axs[location].plot(x, MSError[:,0], '-', color ='black', linewidth='1.5', label = 'Mean error');
    axs[location].plot(x, MSError[:,1], '--', color ='green', linewidth='2', label = '5-Percentile error');
    axs[location].plot(x, MSError[:,2], '-.', color ='red', linewidth='2', label = '95-Percentile error');
    axs[location].plot(x, MSError[:,3], '.', color ='blue', linewidth='3', label = 'Standard deviation');

    axs[location].set_xticks([0,1,2,3,4,5,6]);
    axs[location].set_xticklabels(['12','25','50','100','200','400','800']);

    axs[location].set_ylabel(r"MSE");
    axs[location].legend(loc='upper left', fontsize=5);

    #_______________________________________________________
    # Rel. Min error plot
    #_______________________________________________________
    location = (0,1);
    axs[location].plot(x, MinError[:,0], '-', color ='black', linewidth='1.5', label = 'Mean error');
    axs[location].plot(x, MinError[:,1], '--', color ='green', linewidth='2', label = '5-Percentile error');
    axs[location].plot(x, MinError[:,2], '-.', color ='red', linewidth='2', label = '95-Percentile error');
    axs[location].plot(x, MinError[:,3], '.', color ='blue', linewidth='3', label = 'Standard deviation');

    axs[location].set_xticks([0,1,2,3,4,5,6]);
    axs[location].set_xticklabels(['12','25','50','100','200','400','800']);

    axs[location].set_ylabel(r"Rel. Min Error");
    axs[location].legend(loc='upper left', fontsize=5);

    #_______________________________________________________
    # Max error plot
    #_______________________________________________________
    location = (1,1);
    axs[location].plot(x, MaxError[:,0], '-', color ='black', linewidth='1.5', label = 'Mean error');
    axs[location].plot(x, MaxError[:,1], '--', color ='green', linewidth='2', label = '5-Percentile error');
    axs[location].plot(x, MaxError[:,2], '-.', color ='red', linewidth='2', label = '95-Percentile error');
    axs[location].plot(x, MaxError[:,3], '.', color ='blue', linewidth='3', label = 'Standard deviation');

    axs[location].set_xticks([0,1,2,3,4,5,6]);
    axs[location].set_xticklabels(['12','25','50','100','200','400','800']);

    axs[location].set_ylabel(r"Rel. Max Error");
    axs[location].legend(loc='upper left', fontsize=5);

    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,00))
   
    os.chdir(projectDir);

    plt.savefig("ErrorVariation.pdf");

    os.chdir(curDir);

    pass;


