import numpy as np
from collections import namedtuple

class MuellerBrownPotential(object):

    """Mueller Brown potential."""

    dimension = 2

    # named tuple to store the parameters of the potential
    ParameterTuple = namedtuple('ParameterTuple','A a b c x0 y0')
    param = ParameterTuple(
                    A = np.array( [ -200.0, -100.0, -175.0, 15.0 ] ),
                    a = np.array( [ -1.0, -1.0, -6.5, 0.7 ] ),
                    b = np.array( [ 0.0, 0.0, 11.0, 0.6 ] ),
                    c = np.array( [ -10.0, -10.0, -6.5, 0.7 ] ),
                    x0 = np.array( [ 1.0, 0.0, -0.5, -1.0 ] ),
                    y0 = np.array( [ 0.0, 0.5, 1.5, 1.0 ] )
                )

    minima = (
                np.array( [ -0.558,  1.442 ] ),
                np.array( [  0.623,  0.028 ] ),
                np.array( [ -0.050,  0.467 ] )
            )

    saddlePoints = (
                np.array( [ -0.822,  0.624 ] ),
                np.array( [ -0.212,  0.293 ] )
            )

    gridLimitsMin = [ -1.5, -0.5 ]
    gridLimitsMax = [  1.5,  2.5 ]

    defaultMCStepSize = 0.2

    def __init__(self,scalePotential=0.2):
        # shift the potential such that its value at the global
        # minimum (-0.558,1.442) is 0.0
        self.scalePotentialFactor = 1.0
        self.potentialShift = 0.0
        self.potentialShift = - self.calcMBpotential(self.minima[0])
        self.scalePotentialFactor = scalePotential
        self.potentialPlotLimits = [ 0, 300*self.scalePotentialFactor ]
    #-------------------------------

    def getPotentialValue(self,position):
        potential = self.calcMBpotential(position)
        return potential
    #-------------------------------

    def calcMBpotential(self,position):
        assert len(position) == self.dimension
        value = 0.0
        x = position[0]; y = position[1]
        for i in range(self.param.A.size):
            A = self.param.A[i]
            a = self.param.a[i]; b = self.param.b[i]; c = self.param.c[i]
            x0 = self.param.x0[i]; y0 = self.param.y0[i]
            tmp1 = a*(x-x0)**2 + b*(x-x0)*(y-y0) + c*(y-y0)**2
            value += A*np.exp(tmp1)
        value += self.potentialShift
        value = value * self.scalePotentialFactor
        return value
    #-------------------------------

    def getMinima(self):
        return self.minima
    #-------------------------------

    def getSaddlePoints(self):
        return self.saddlePoints
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
