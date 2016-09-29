import time
import numpy as np
from scipy import special
import matplotlib.pyplot as plt
pi = np.pi

def ExactValueOfSphereVolume(Dimension = 3, Radius = 1.0):
    Dim = Dimension
    ExactValue = ( pi**(Dim/2.0) * Radius**Dim ) / special.gamma( Dim/2.0+1.0 )
    return ExactValue
#-----------------------

def getEstimateOfSphereVolume(NumSamples, Dimension = 3, Radius = 1.0, InitialSeed = -1):
    Dim = Dimension
    if InitialSeed >= 0:
        np.random.seed(InitialSeed)
    else:
        np.random.seed()
    BoxVolume = (2.0 * Radius)**Dim
    estimate = []
    inside = 0
    total = 0
    for i in range(0,NumSamples):
        sum = 0.0
        for d in range(Dim):
            x = np.random.uniform(-Radius,Radius)
            sum += x*x
        sum = np.sqrt(sum)
        if(sum <= Radius): inside += 1
        total += 1
        ratio = np.float(inside)/np.float(total)
        value = ratio * BoxVolume
        estimate.append(value)
    estimate = np.array(estimate)
    return np.copy(estimate)
#----------------------------


NumRuns = 20
NumSamples = 1000000
Dimension = 5
ExactValue = ExactValueOfSphereVolume(Dimension = Dimension)

FigOffset = 0.5

y_min = ExactValue - FigOffset
y_max = ExactValue + FigOffset
plt.axis([0, NumSamples, y_min, y_max])

for i in range(NumRuns):
    estimate = getEstimateOfSphereVolume(NumSamples = NumSamples, Dimension = Dimension)
    plt.plot(estimate)
    if (i+1) % 10 == 0: print 'Run {:4d} of {:4d}'.format(i+1,NumRuns)
plt.plot([0, NumSamples], [ExactValue, ExactValue], 'r--')
plt.show()
