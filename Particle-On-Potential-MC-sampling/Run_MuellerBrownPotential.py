import numpy as np
import time
import matplotlib.pyplot as plt
from MuellerBrownPotential import MuellerBrownPotential
from LogExpOfHarmonicWellsPotential import LogExpOfHarmonicWellsPotential
from MonteCarloSimulator import MonteCarloSimulator

T = 1.0
NumMCmoves = 10000

# potential = LogExpOfHarmonicWellsPotential()
potential = MuellerBrownPotential()
MCsim = MonteCarloSimulator(
            potentialClass=potential,
            Temperature = T
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
