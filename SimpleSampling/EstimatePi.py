import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

pi = np.pi

def getEstimateOfPi(NumSamples,InitialSeed=-1):
    if InitialSeed >= 0:
        np.random.seed(InitialSeed)
    else:
        np.random.seed()
    estimates=[]
    inside = 0
    total = 0
    for i in range(0,NumSamples):
        x = np.random.uniform(-1.0,1.0)
        y = np.random.uniform(-1.0,1.0)
        if x**2+y**2 <= 1.0**2:
            inside += 1
        total += 1
        current_estimate = 4 * np.float(inside)/np.float(total)
        estimates.append( current_estimate )
    estimates = np.array(estimates)
    return np.copy(estimates)
#----------------------------------

def getEstimateOfPiIntegral(NumSamples,InitialSeed=-1):
    if InitialSeed >= 0:
        np.random.seed(InitialSeed)
    else:
        np.random.seed()
    estimates=[]
    current_estimate = 0.0
    for i in range(0,NumSamples):
        x = np.random.uniform(0,1.0)
        value = 4.0 * np.sqrt(1-x**2)
        current_estimate += (value - current_estimate)/(float(i)+1.0)
        estimates.append( current_estimate )
    estimates = np.array(estimates)
    return np.copy(estimates)
#----------------------------------

def getEstimateOfPiMarkovChain(NumSamples,InitialSeed=-1,MaxStepSize=0.1):
    if InitialSeed >= 0:
        np.random.seed(InitialSeed)
    else:
        np.random.seed()
    estimates=[]
    inside = 0
    total = 0
    accepted = 0
    x = np.random.uniform(-1.0,1.0)
    y = np.random.uniform(-1.0,1.0)
    for i in range(0,NumSamples):
        dx = np.random.uniform(-MaxStepSize,MaxStepSize)
        dy = np.random.uniform(-MaxStepSize,MaxStepSize)
        x_new = x + dx
        y_new = x + dx
        if abs(x_new) < 1.0 and abs(y_new) < 1.0:
            x = x_new
            y = y_new
            accepted += 1
        if x**2+y**2 <= 1.0**2:
            inside += 1
        total += 1
        current_estimate = 4 * np.float(inside)/np.float(total)
        estimates.append( current_estimate )
    estimates = np.array(estimates)
    acceptance_ratio = np.float(accepted)/np.float(total)
    return np.copy(estimates), acceptance_ratio
#----------------------------------

def runManyRuns(NumSamples,NumRuns=50,function=getEstimateOfPi):
    AllEstimates = []
    for i in range(NumRuns):
        estimate = function(NumSamples)
        AllEstimates.append(estimate)
    AllEstimates = np.array(AllEstimates)
    N = np.arange(1,NumSamples+1)
    SqrtN = np.sqrt(N)
    MeanN = []
    StddevN = []
    for j in range(NumSamples):
        mean = np.mean(AllEstimates[:,j])
        stddev = np.std(AllEstimates[:,j])
        MeanN.append(mean)
        StddevN.append(stddev)
    MeanN = np.array(MeanN)
    StddevN = np.array(StddevN)
    Data = (
        np.copy(N),
        np.copy(SqrtN),
        np.copy(MeanN),
        np.copy(StddevN),
        np.copy(AllEstimates)
    )
    return Data
#----------------------------------

def plotData(Data, FigNumber=1):
    exact_value = pi
    #
    N = Data[0]
    SqrtN = Data[1]
    MeanN = Data[2]
    StddevN = Data[3]
    AllEstimates = Data[4]
    NumSamples = AllEstimates[0,:].size
    NumRuns = AllEstimates[:,0].size
    #
    FigOffset1 = 0.3
    plt.figure(FigNumber)
    plt.xlim(0, NumSamples)
    plt.ylim(exact_value-FigOffset1, exact_value+FigOffset1)
    plt.xticks( np.linspace(0,NumSamples,6, endpoint=True) )
    # plt.yticks(
    #    [pi-1.2,pi-1.0, pi-0.8, pi-0.6, pi-0.4, pi-0.2, pi-0.0, pi+0.2, pi+0.4, pi+0.6, pi+0.8, pi+1.0, pi+1.2],
    #    [r'$\pi-1.2$', r'', r'$\pi-0.8$', r'', r'$\pi-0.4$', r'', r'$\pi$', r'',
    #     r'$\pi+0.4$', r'', r'$\pi+0.8$', r'', r'$\pi+1.2$'],
    #    fontsize = 20)
    for i in range(NumRuns):
        plt.plot( AllEstimates[i,:] )
    plt.plot([0, NumSamples], [exact_value, exact_value], 'k--', linewidth=4)
    #
    plt.figure(FigNumber+1)
    plt.plot(N,StddevN, 'b--', linewidth=4)
    plt.figure(FigNumber+2)
    plt.plot(N,SqrtN*StddevN, 'b--', linewidth=4)
    plt.show()
#----------------------------------




#(N, sqrtN, MeanN, StddevN, AllEstimates) = runManyRuns(10000, NumRuns=100, function=getEstimateOfPiIntegral)
(N, sqrtN, MeanN, StddevN, AllEstimates) = runManyRuns(10000, NumRuns=1000, function=getEstimateOfPi)
Data = (N, sqrtN, MeanN, StddevN, AllEstimates)
plotData(Data)








# mean = np.mean(FinalValues)
# stddev = np.std(FinalValues)
# min = np.min(FinalValues)
# max = np.max(FinalValues)
#
#
# plt.figure(2)
# plt.hist(FinalValues,bins=40,normed=True,histtype='bar')
# x = np.linspace(min,max,100)
# plt.plot(x, norm.pdf(x, loc=mean, scale=stddev))
# print x
# print norm.pdf(x, loc=mean, scale=stddev)
# plt.show()
