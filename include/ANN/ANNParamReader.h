// =======================================================================
// @(#)TANNParamReader.h        1.0 24.12.20
// 
// Class:       TANNParamReader
// Purpose:     Class for reading the parameter file 
//
// Author:      Subodh M. Joshi (28.12.20)
//
// History:     start of implementation 28.12.20 (Subodh M. Joshi)
//
// =======================================================================

#ifndef __TANNPARAMREADER__
#define __TANNPARAMREADER__

#include <ANNFunctions.h>

#ifdef _MPI
#include "mpi.h"
#endif

class TANNParamReader
{
  public:
    /** constructor */
    TANNParamReader(char *paramFile);

    /** destructor */
    ~TANNParamReader();

    /** number of hidden layers in FFN */
    int nHL;

    /** number of inputs to each input neuron. i.e. size of the data per perceptron at input. This could be considered as r,g,b,intensity values per neuron at input. */
    int ipDataDim;

    /** array storing dim of each layer (from IP to OP) */
    int *layerDim;

    /** Array storing types of each layer i.e. activation functions
     * NOTE: for the input layer, the activation function is 'dummy'
     * */
    std::string *layerType;

    /** array storing int code for the layer type */
    int *layerTypeInt;

    /** Int code storing the error type (loss function) */
    int lossFunctionInt;

    /** loss function */
    std::string lossFunction;

    /** Training Dataset name */
    std::string trainingDatasetName;

    /** Testing Dataset name */
    std::string testingDatasetName;

    /** Validation data percentage of the total training data */
    double validationDataPercentage;

    /** Optimizer code */
    /** 0: Default (RMSProp)
     *  1: Gradient Descent
     *  2: SGD
     *  3: Adam
     *  4: L-BFGS
     *  **/
    int optimizerCode;

    /** Step size for the optimizer */
    double optimizerStepSize;

    /** Stochastic GD batch size */
    int sgdBatchSize;

    /** Max number of iterations */
    int maxIterations;

    /** Constant used for feature scaling */
    double featureScalingConstant;

    int epochs;

    // Parameter used for regularization before the output layer
    double dropoutRatio;

    // Tolerance for the optimizer
    double tolerance;

    // File for saving the model
    std::string saveModelFile;

    // File for loading the model
    std::string loadModelFile;

    // Loading vs training flag
    // This is set to 1 if the model is read from a file instead of training from scratch
    int loadModelFlag;

    // File for saving the data
    std::string saveDataFile;

  private:
    /** flags */
    bool layerDimFlag;
    bool layerTypeFlag;
    bool layerFlag;
    
  public:
    /** Class methods **/

    /** read total number of layers */
    void readNumberOfLayers(char *paramFile);

    /** read param file */
    void readParamFile(char *paramFile);

    /** Print the data */
    void print();

  private:
    /** Class methods private */
    std::string getLayerType(int number);

    std::string getLossFunction(int number);
};

#endif
