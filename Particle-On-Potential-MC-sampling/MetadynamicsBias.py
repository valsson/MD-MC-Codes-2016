import numpy as np
import matplotlib.pyplot as plt


class MetadynamicsBias(object):

    """
       Code to add a metadynamics bias in N-dimensions.
    """

    pi = np.pi
    kB = 1.0

    gridBinsDefault = 200

    def __init__(self,
                Temperature,
                Sigma,
                InitialHeight,
                Pace=500,
                Biasfactor=-1
                ):
        self.T = float(Temperature)
        self.kBT = self.kB * self.T
        self.beta = 1.0/self.kBT
        # Metadynamics parameters
        self.sigma = np.array(Sigma)
        self.dimension = len(self.sigma)
        self.initialHeight = InitialHeight
        self.pace = Pace
        self.wellTempered = False
        self.gamma = Biasfactor
        if self.gamma >= 1.0:
            self.wellTempered = True
        #
        self.hillCenters = []
        self.hillHeights = []
        self.hillCenters.append( np.zeros(self.dimension) )
        self.hillHeights.append( 0.0 )
    #-------------------------------

    def getBiasValue(self,cv):
        # assert len(cv) == self.dimension
        centersTmp = np.array(self.hillCenters)
        heightsTmp = np.array(self.hillHeights)
        tmpArg = np.zeros(heightsTmp.shape)
        for k in xrange(self.dimension):
            tmpArg += - np.power( cv[k]-centersTmp[:,k], 2 ) / ( 2.0*self.sigma[k]**2 )
        gaussiansTmp = heightsTmp * np.exp(tmpArg)
        bias = np.sum(gaussiansTmp)
        #bias = 0.0
        return bias
    #-------------------------------

    def addGaussianHill(self,cv):
        # assert len(cv) == self.dimension
        currentCenter = np.array(cv)
        sigma = self.sigma
        currentBias = self.getBiasValue(currentCenter)
        WTscalingFactor = 1.0
        if self.wellTempered:
            WTscalingFactor = np.exp(- (1.0/(self.gamma-1.0)) * self.beta * currentBias )
        currentHeight = self.initialHeight * WTscalingFactor
        self.hillCenters.append( currentCenter )
        self.hillHeights.append( currentHeight )
    #-------------------------------

    def updateBias(self,step,cv):
        if step % self.pace == 0:
            self.addGaussianHill(cv)
    #-------------------------------

    def getGaussianCenters(self):
        return np.array(self.hillCenters)
    #-------------------------------

    def getGaussianHeights(self):
        return np.array(self.hillHeights)
    #-------------------------------

    def getCentersMinMax(self):
        centers = self.getGaussianCenters()
        centersMin = np.zeros(self.dimension)
        centersMax = np.zeros(self.dimension)
        for i in range(self.dimension):
            centersMin[i] = np.min(centers[:,i])
            centersMax[i] = np.max(centers[:,i])
        return (centersMin, centersMax)
    #-------------------------------

    def getEstimatedFES(self,gridBins=[],gridMin=[],gridMax=[]):
        if len(gridBins) == 0:
            gridBins = np.full(self.dimension, self.gridBinsDefault, dtype=np.int64)
        assert len(gridBins)==self.dimension
        if len(gridMin) and len(gridMax) == 0:
            (gridMin,gridMax) = self.getCentersMinMax(self)
        assert len(gridMin) == self.dimension
        assert len(gridMax) == self.dimension
        if self.dimension == 1:
            (xGrid, fesGrid) = self.calculateFES_1D(gridBins,gridMin,gridMax)
            return (xGrid, fesGrid)
        elif self.dimension == 2:
            (xGrid, yGrid, fesGrid) = self.calculateFES_2D(gridBins,gridMin,gridMax)
            return (xGrid, yGrid, fesGrid)
        else:
            pass
    #-------------------------------

    def calculateFES_1D(self,gridBins,gridMin,gridMax,minToZero=True):
        assert len(gridBins) == self.dimension
        assert len(gridMin) == self.dimension
        assert len(gridBins) == self.dimension
        xGrid = np.linspace(gridMin[0], gridMax[0], gridBins[0])
        fesGrid = np.zeros(gridBins)
        wtFactor = 1.0
        if self.wellTempered:
            wtFactor = self.gamma/(self.gamma-1.0)
        for i in range(gridBins[0]):
            fesGrid[i] = - wtFactor * self.getBiasValue([ xGrid[i] ])
        if(minToZero):
            fesGrid = fesGrid - np.min(fesGrid)
        return (xGrid, fesGrid)
    #-------------------------------

    def calculateFES_2D(self,gridBins,gridMin,gridMax,minToZero=True):
        assert len(gridBins) == self.dimension
        assert len(gridMin) == self.dimension
        assert len(gridMax) == self.dimension
        xTmpGrid = np.linspace(gridMin[0], gridMax[0], gridBins[0])
        yTmpGrid = np.linspace(gridMin[1], gridMax[1], gridBins[1])
        [ xGrid, yGrid ] = np.meshgrid(xTmpGrid,yTmpGrid)
        fesGrid = np.zeros(gridBins)
        wtFactor = 1.0
        if self.wellTempered:
            wtFactor = self.gamma/(self.gamma-1.0)
        for i in range(gridBins[0]):
            for j in range(gridBins[1]):
                fesGrid[i,j] = - wtFactor * self.getBiasValue([ xGrid[i,j], yGrid[i,j] ])
        if(minToZero):
            fesGrid = fesGrid - np.min(fesGrid)
        return (xGrid, yGrid, fesGrid)
    #-------------------------------

    def plotFES_1D(self,gridBins,gridMin,gridMax,yLim,FigNumber,FigFilename=None):
        plt.figure(FigNumber)
        plt.clf()
        (xGrid, fesGrid) = self.calculateFES_1D(gridBins,gridMin,gridMax)
        plt.plot(xGrid,fesGrid,lw=3, color='black')
        plt.ylim(yLim)
        plt.xlabel('x')
        plt.ylabel('F(x)')
        if FigFilename != None:
            FigFormat = FigFilename[-3:]
            plt.savefig(FigFilename, transparent=True, format=FigFormat)
        plt.show()
    #-------------------------------

    def plotFES_2D(self,gridBins,gridMin,gridMax,zLim,FigNumber,FigFilename=None):
        plt.figure(FigNumber)
        plt.clf()
        (xGrid, yGrid, fesGrid) = self.calculateFES_2D(gridBins,gridMin,gridMax)
        plt.contourf(xGrid, yGrid, fesGrid,
                            vmin=zLim[0], vmax=zLim[1],
                            levels=np.linspace(zLim[0],zLim[1],1000) )
        plt.colorbar(ticks=np.linspace(zLim[0],zLim[1],5), boundaries=np.linspace(-zLim[0],zLim[1],5))
        plt.xlabel('x')
        plt.ylabel('y')
        if FigFilename != None:
            FigFormat = FigFilename[-3:]
            plt.savefig(FigFilename, transparent=True, format=FigFormat)
        plt.show()
    #-------------------------------

    def getTemperature(self):
        return self.T
    #-------------------------------

    def setTemperature(self,value):
        self.T = value
    #-------------------------------

    def getBeta(self):
        return self.beta
    #-------------------------------
