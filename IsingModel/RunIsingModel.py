import numpy as np
import time
from IsingModel import IsingModel
import matplotlib.pyplot as plt
from DataTools import *

NumSweeps = 10000

plotSpinStateEvery = 10
m1 = IsingModel(FigNum=1, L=11, Temperature=3.0, Dimension=2,InitialState='random')
# m1.randomizeSpins()

energyPerSpin = []
magnetizationPerSpin = []
absMagnetizationPerSpin = []

energyPerSpin.append( m1.getEnergyPerSpin() )
magnetizationPerSpin.append( m1.getMagnetizationPerSpin() )
absMagnetizationPerSpin.append( m1.getAbsoluteMagnetizationPerSpin() )

for i in range(NumSweeps):
    m1.doOneMonteCarloSweep()
    energyPerSpin.append( m1.getEnergyPerSpin() )
    magnetizationPerSpin.append( m1.getMagnetizationPerSpin() )
    absMagnetizationPerSpin.append( m1.getAbsoluteMagnetizationPerSpin() )
    if (i+1) % plotSpinStateEvery == 0:
        m1.updateSpinFig()

energyPerSpin = np.array(energyPerSpin)
magnetizationPerSpin = np.array(magnetizationPerSpin)
absMagnetizationPerSpin = np.array(absMagnetizationPerSpin)

data = [energyPerSpin,magnetizationPerSpin,absMagnetizationPerSpin]
dataFields = ['energyPerSpin','magnetizationPerSpin','absMagnetizationPerSpin']

# writeDataToFile('Ising.T30.data', data, fieldNames=dataFields, addTimeField=True)

plt.figure(2)
plt.plot(energyPerSpin)
plt.draw()
plt.ylim(-2.1,0.6)
plt.figure(3)
plt.plot(magnetizationPerSpin)
plt.ylim(-1.1,1.1)
plt.draw()
plt.figure(4)
plt.plot(absMagnetizationPerSpin)
plt.ylim(-0.1,1.1)
plt.draw()
time.sleep(100.0)
