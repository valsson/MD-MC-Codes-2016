import numpy as np
import time
import matplotlib.pyplot as plt
from LogExpOfHarmonicWellsPotential import LogExpOfHarmonicWellsPotential
from MonteCarloSimulator import MonteCarloSimulator

MBpot = LogExpOfHarmonicWellsPotential()
MCsim = MonteCarloSimulator(
                        potentialClass=MBpot,
                        Temperature = 1.0
                        )

StepSizes = np.linspace(5.1,10.0,50)

Data = []
for i in range(StepSizes.size):
    MCsim.resetRun()
    StepSize = StepSizes[i]
    MCsim.setStepSize( StepSize )
    MCsim.setPosition( MBpot.getMinima()[0] )
    MCsim.doMonteCarloMoves(100000,verbose=False)
    # MCsim.plotTrajectoryOnPotential()
    AveAcceptence = MCsim.getAverageAcceptence()
    (Mean, StdDev) = MCsim.getTrajectoryMeanAndStdDev()
    AveDisplacement = MCsim.getAverageDisplacement()
    Data.append( [StepSize, AveAcceptence, Mean, StdDev[0], AveDisplacement] )
Data = np.array(Data)

plt.figure(2)
plt.plot(
        Data[:,0],
        Data[:,1],
        lw=4,
        marker='x',
        markersize=20,
        markeredgecolor='red'
        )

# plt.savefig('StepSizes.pdf',format='pdf',transparent=True)

plt.figure(3)
plt.plot(
        Data[:,0],
        Data[:,2],
        lw=4,
        marker='x',
        markersize=20,
        markeredgecolor='red'
        )

plt.figure(4)
plt.plot(
        Data[:,0],
        Data[:,3],
        lw=4,
        marker='x',
        markersize=20,
        markeredgecolor='red'
        )

plt.figure(5)
plt.plot(
        Data[:,0],
        Data[:,4],
        lw=4,
        marker='x',
        markersize=20,
        markeredgecolor='red'
        )

plt.show()
