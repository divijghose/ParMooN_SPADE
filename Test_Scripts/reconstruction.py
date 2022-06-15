from fileinput import filename
import numpy as np 
import pandas as pd 
from matplotlib import pyplot as  plt 
import seaborn as sns
import os
import pickle

InData = pd.read_csv("PyIn.txt",',',header=None)
N_R = InData[0][0]
subDim = InData[0][1]
maxTimeStep = InData[0][2]

meandiffmax=[]
meandiffmin = []
meandiffmean = []
meandiffmse = []
meandiffrmse = []

recondiffmax=[]
recondiffmin = []
recondiffmean = []
recondiffmse = []
recondiffrmse = []

recondiffrmsemax = []
recondiffrmsemin = []
recondiffrmsemean = []


currdir = os.getcwd()
errdir = os.path.join(currdir,"Error_Analysis")
if(not os.path.exists(errdir)):
	os.mkdir(errdir)
    
currdir = os.getcwd()
pkldir = os.path.join(currdir,"Error_Analysis/Error_Data")
if(not os.path.exists(pkldir)):
	os.mkdir(pkldir)
    
currdir = os.getcwd()
eplotdir = os.path.join(currdir,"Error_Analysis/Error_Plots")
if(not os.path.exists(eplotdir)):
	os.mkdir(eplotdir)

for i in range(maxTimeStep):
    if(i<10):
        filenameCoeff = "Coefficients/Coeff_NRealisations_"+str(N_R)+"_t0000" + str(i)+".txt"
    elif(i<100):
        filenameCoeff = "Coefficients/Coeff_NRealisations_"+str(N_R)+"_t000" + str(i)+".txt"
    elif(i<1000):
        filenameCoeff = "Coefficients/Coeff_NRealisations_"+str(N_R)+"_t00" + str(i)+".txt"
    else:
        filenameCoeff = "Coefficients/Coeff_NRealisations_"+str(N_R)+"_t0" + str(i)+".txt"
        
    if(i<10):
        filenameMean = "Mean/Mean_NRealisations_"+str(N_R)+"_t0000" + str(i)+".txt"
    elif(i<100):
        filenameMean = "Mean/Mean_NRealisations_"+str(N_R)+"_t000" + str(i)+".txt"
    elif(i<1000):
        filenameMean = "Mean/Mean_NRealisations_"+str(N_R)+"_t00" + str(i)+".txt"
    else:
        filenameMean = "Mean/Mean_NRealisations_"+str(N_R)+"_t0" + str(i)+".txt"
        
    if(i<10):
        filenameMode = "Modes/Mode_NRealisations_"+str(N_R)+"_t0000" + str(i)+".txt"
    elif(i<100):
        filenameMode = "Modes/Mode_NRealisations_"+str(N_R)+"_t000" + str(i)+".txt"
    elif(i<1000):
        filenameMode = "Modes/Mode_NRealisations_"+str(N_R)+"_t00" + str(i)+".txt"
    else:
        filenameMode = "Modes/Mode_NRealisations_"+str(N_R)+"_t0" + str(i)+".txt"
        
    if(i<10):
        filenameMC = "MonteCarlo/MC_NRealisations_"+str(N_R)+"_t0000" + str(i)+".txt"
    elif(i<100):
        filenameMC = "MonteCarlo/MC_NRealisations_"+str(N_R)+"_t000" + str(i)+".txt"
    elif(i<1000):
        filenameMC = "MonteCarlo/MC_NRealisations_"+str(N_R)+"_t00" + str(i)+".txt"
    else:
        filenameMC = "MonteCarlo/MC_NRealisations_"+str(N_R)+"_t0" + str(i)+".txt"
        
    coeffarr = pd.read_csv(filenameCoeff,',',header=None)
    modearr = pd.read_csv(filenameMode,',',header=None)
    meanarr = pd.read_csv(filenameMean,header=None)
    mcarr = pd.read_csv(filenameMC,',',header=None)
    
    A = mcarr.to_numpy() #Monte Carlo 
    B = (meanarr[0].to_numpy()[:,np.newaxis]+np.matmul(modearr.to_numpy(),np.transpose(coeffarr.to_numpy()))) #Reconstructed Matrix
    C = A.mean(axis=1)[:,np.newaxis] #Monte Carlo Mean
    D = meanarr[0].to_numpy()[:,np.newaxis] #Mean solution
    meandiff = np.abs(C-D) #Absolute difference between Monte Carlo mean and mean solution
    meandiffmax.append(meandiff.max())
    meandiffmin.append(meandiff.min())    
    meandiffmean.append(meandiff.mean())
    
    meandiffmse.append((np.square(C - D)).mean(axis=0))#mean squared error between MC mean and mean solution
    meandiffrmse.append(np.sqrt((np.square(C - D)).mean(axis=0)))#root mean squared error between MC mean and mean solution
    
    recondiff = np.abs(A-B) #Absolute difference between Monte Carlo mean and mean solution
    
    recondiffmax.append(recondiff.max()) #max absolute error
    filename = 'Error_Analysis/Error_Data/recondiffmax_NR'+str(N_R)+'.pkl'
    with open(filename, 'wb') as f:
        pickle.dump(recondiffmax,f)
        
#     recondiffmin.append(recondiff.min())  
    
    recondiffmean.append(recondiff.mean())
    filename = 'Error_Analysis/Error_Data/recondiffmean_NR'+str(N_R)+'.pkl'
    with open(filename, 'wb') as f:
        pickle.dump(recondiffmean,f)
    
    recondiffmse.append(np.max((np.square(A - B)).mean(axis=0)))#mean squared error between MC mean and mean solution
    filename = 'Error_Analysis/Error_Data/recondiffmse_NR'+str(N_R)+'.pkl'
    with open(filename, 'wb') as f:
        pickle.dump(recondiffmse,f)
        
    recondiffrmsemax.append(np.max(np.sqrt((np.square(A - B)).mean(axis=0))))#max of root mean squared error between MC mean and mean solution
    filename = 'Error_Analysis/Error_Data/recondiffrmsemax_NR'+str(N_R)+'.pkl'
    with open(filename, 'wb') as f:
        pickle.dump(recondiffrmsemax,f)
        
    recondiffrmsemin.append(np.min(np.sqrt((np.square(A - B)).mean(axis=0))))#max of root mean squared error between MC mean and mean solution
    filename = 'Error_Analysis/Error_Data/recondiffrmsemin_NR'+str(N_R)+'.pkl'
    with open(filename, 'wb') as f:
        pickle.dump(recondiffrmsemin,f)
    recondiffrmsemean.append(np.mean(np.sqrt((np.square(A - B)).mean(axis=0))))#max of root mean squared error between MC mean and mean solution
    filename = 'Error_Analysis/Error_Data/recondiffrmsemean_NR'+str(N_R)+'.pkl'
    with open(filename, 'wb') as f:
        pickle.dump(recondiffrmsemean,f)
    
    
    
    
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['axes.labelsize']=15

plt.rcParams["font.serif"]
SMALL_SIZE = 15
MEDIUM_SIZE = 15
BIGGER_SIZE = 15

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title



plt.plot(meandiffmax, 'r--',label="Maximum absolute error")
# p = plt.plot(meandiffmin, 'b-',label="Minimum absolute error")
plt.plot(meandiffmean, 'g-',label="Mean absolute error")
plt.plot(meandiffrmse, 'b-x',label="Root mean squared error")

plt.xlabel("Time-step")
plt.ylabel("Error")
plt.legend()

plt.title("Error between solution of mean equation\nand mean of Monte Carlo solutions, No. of Realizations: "+str(N_R))
plt.gcf().set_size_inches(11, 8)

plt.savefig("Error_Analysis/Error_Plots/MeanError_NR"+str(N_R)+".png", facecolor='white',bbox_inches='tight')
plt.close()

plt.rcParams['xtick.labelsize'] = 24
plt.rcParams['ytick.labelsize'] = 24
plt.rcParams['axes.labelsize']=24

plt.rcParams["font.serif"]
SMALL_SIZE = 24
MEDIUM_SIZE = 24
BIGGER_SIZE = 24

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


fig,ax=plt.subplots(2,2,sharex=True,sharey=False)

ax[0,0].plot(recondiffrmsemax, 'b-x',label="Maximum RMSE")
ax[0,0].plot(recondiffrmsemin, 'g--',label="Minimum RMSE")
ax[0,0].plot(recondiffrmsemean, 'r-',label="Mean RMSE")
ax[0,0].set_ylabel("Error\n")
# ax[0,0].legend()
ax[0,1].plot(recondiffrmsemax, 'b-x')

ax[1,0].plot(recondiffrmsemin, 'g--')
ax[1,0].set_xlabel("Time-step")
ax[1,0].set_ylabel("Error")



ax[1,1].plot(recondiffrmsemean, 'r-')
ax[1,1].set_xlabel("Time-step")
ax[0,0].grid()
ax[0,1].grid()
ax[1,0].grid()
ax[1,1].grid()
fig.legend(ncol=3,loc=9,bbox_to_anchor=(0.5,0.95))
fig.suptitle("Reconstruction error across realizations, No. of Realizations: "+str(N_R))
plt.gcf().set_size_inches(22, 16)
plt.savefig("Error_Analysis/Error_Plots/Recon_RMSE_NR"+str(N_R)+".png", facecolor='white',bbox_inches='tight')
plt.close()
