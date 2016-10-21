import numpy as np
import time
import matplotlib.pyplot as plt
from MuellerBrownPotential import MuellerBrownPotential
from LogExpOfHarmonicWellsPotential import LogExpOfHarmonicWellsPotential
from MonteCarloSimulator import MonteCarloSimulator
from MetadynamicsBias import MetadynamicsBias


T = 1.0
NumMCmoves = 10000

kB = 1.0
kBT = kB * T

initialHeight = (kBT/2.0)
sigma = [ 0.2,0.2 ]
pace = 1
biasfactor = 10.0

# potential = LogExpOfHarmonicWellsPotential()
potential = MuellerBrownPotential()
MetaD = MetadynamicsBias(
            Temperature = T,
            Sigma = sigma,
            InitialHeight = initialHeight,
            Pace = pace,
            Biasfactor = biasfactor
            )
MCsim = MonteCarloSimulator(
            potentialClass=potential,
            Temperature = T,
            externalBiasClass = MetaD
            )


MCsim.resetRun()
MCsim.setPosition( potential.getMinima()[0] )
MCsim.runMC(NumMCmoves)
print ' '
MCsim.printAverageAcceptence()
print ' '
MCsim.printTrajectoryMeanAndStddev()
print ' '
MCsim.plotPotentialAndTrajectory()
MCsim.plotTrajectoryTimeSeries()
MCsim.plotTrajectoryHistogramAndFES()
