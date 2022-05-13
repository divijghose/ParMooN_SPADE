from fileinput import filename
import numpy as np 
import pandas as pd 
from matplotlib import pyplot as  plt 
import seaborn as sns

InData = pd.read_csv("PyIn.txt",',',header=None)
N_R = InData[0][0]
subDim = InData[0][1]
maxTimeStep = InData[0][2]

coeffnamearr =[]
for i in range(subDim):
    coeffname = "Coefficient "+str(i+1)
    coeffnamearr.append(coeffname)

for i in range(maxTimeStep):
    filename = "Coefficients/Coeff_NRealisations_"+str(N_R)+"_t" + str(i)
    coeffarr = pd.read_csv(filename,',',header=None,names=coeffnamearr)

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
    p=sns.pairplot(data=coeffarr,kind="kde")
    p.fig.suptitle("Joint Probability Distribution, Timestep : "+str(i))
    p.fig.subplots_adjust(top=0.9)
    p.fig.set_size_inches(11,8)
    plt.savefig("PairPlot_t"+str(i)+".png")

