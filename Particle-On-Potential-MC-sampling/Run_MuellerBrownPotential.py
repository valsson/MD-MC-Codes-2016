import numpy as np
import time
import matplotlib.pyplot as plt
from MuellerBrownPotential import MuellerBrownPotential
from LogExpOfHarmonicWellsPotential import LogExpOfHarmonicWellsPotential
from MonteCarloSimulator import MonteCarloSimulator

T = 2.0

# potential = LogExpOfHarmonicWellsPotential()
potential = MuellerBrownPotential()
MCsim = MonteCarloSimulator(
            potentialClass=potential,
            Temperature = T
            )


MCsim.resetRun()
MCsim.setPosition( potential.getMinima()[0] )
MCsim.runMC(100000)
MCsim.plotPotentialAndTrajectory(FigFilename='T4.png')
MCsim.plotTrajectoryTimeSeries()
MCsim.plotTrajectoryHistogramAndFES()
MCsim.printAverageAcceptence()
MCsim.printTrajectoryMeanAndStddev()
