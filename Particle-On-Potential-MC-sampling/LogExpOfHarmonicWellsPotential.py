import numpy as np
from collections import namedtuple

class LogExpOfHarmonicWellsPotential(object):

    """
        One-dimensional potential composed of two harmonic potentials that
        are joined by taking the logarithm of a sum of two expoentials that
        correponds to each of the harmonic potentials.

        V(x) = -(1/s)*log( exp(-s*V_a(x)) + exp(-s*V_b(x)) )
        where s is a smoothness parameter and
        V_a(x) = k_a*(x-x0_a) + s_a
        V_b(x) = k_b*(x-x0_b) + s_b
    """

    dimension = 1

    # named tuple to store the parameters of the potential
    ParameterTuple = namedtuple('ParameterTuple','center kappa shift')
    param = ParameterTuple(
                    center = np.array( [ -50.0, +50.0 ] ),
                    kappa = np.array( [  0.05,  0.05 ] ),
                    shift = np.array( [ 0.0, +0.0] )
                )
    gridLimitsMin = [ -100.0 ]
    gridLimitsMax = [  150.0 ]

    defaultMCStepSize = 10.0

    def __init__(self,smoothness=0.01):
        self.potentialPlotLimits = [ -20, 200]
        self.smoothness = float(smoothness)
    #-------------------------------

    def getPotentialValue(self,position):
        potential = self.calcPotential(position)
        return potential
    #-------------------------------

    def calcPotential(self,position):
        assert len(position) == self.dimension
        value = 0.0
        x = position[0]
        smth = self.smoothness
        for i in range(self.param.center.size):
            x0 = self.param.center[i]
            k = self.param.kappa[i]
            shift = self.param.shift[i]
            value += np.exp( - smth * (k * (x-x0)**2 + shift) )
        # value += self.potentialShift
        value = -(1.0/smth)* np.log(value)
        return value
    #-------------------------------

    def getMinima(self):
        minima = ()
        for c in self.param.center:
            minima = minima + (np.array([c]),)
        return minima
    #-------------------------------

    def getDimension(self):
        return self.dimension
    #-------------------------------

    def getGridLimitsMin(self):
        return self.gridLimitsMin
    #-------------------------------

    def getGridLimitsMax(self):
        return self.gridLimitsMax
    #-------------------------------

    def getPlotLimits(self):
        return self.potentialPlotLimits
    #-------------------------------

    def getDefaultMCStepSize(self):
        return self.defaultMCStepSize
    #-------------------------------
