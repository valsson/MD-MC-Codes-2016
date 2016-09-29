import numpy as np
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt

class MonteCarloSimulator(object):

    """
       Code to run Monte Carlo simulations for a one particle on a
       N-dimensional potential. Also has the ability to add a
       metadynamics bias.
    """

    pi = np.pi
    kB = 1.0

    particleMarkerColor   = 'red'
    particleMarker   = 'x'
    particleMarkerSize   = 10.0

    numberOfGridBins = 200

    def __init__(self,
                potentialClass,
                Temperature = 1.0,
                MonteCarloStepSize = None,
                InitialSeed = -1,
                externalBiasClass = None):
        # The potential model used
        self.potentialClass = potentialClass
        self.dimension = self.potentialClass.getDimension()
        # Temperature
        self.T = float(Temperature)
        self.kBT = self.kB * self.T
        self.beta = 1.0/self.kBT
        # Various things
        if InitialSeed >= 0:
            np.random.seed(InitialSeed)
        # Parameters of the Monte Carlo moves
        if MonteCarloStepSize == None:
            self.MCstepSize = potentialClass.getDefaultMCStepSize()
        else:
            self.MCstepSize = float(MonteCarloStepSize)
        # Initialize counters
        self.movesTotal = 0
        self.movesAccepted = 0
        self.movesRejected = 0
        #
        self.totalDisplacement = 0.0
        #
        self.trajectory = []
        self.potentialValues = []
        self.externalBiasValues = []
        # Initialize grid of potential used for plotting
        self.initPotentialGrid()
        #
        self.useExternalBias = False
        self.externalBiasClass = externalBiasClass
        if self.externalBiasClass != None:
            self.useExternalBias = True
        # Initialize position
        tmpPosition = self.randomPosition()
        self.setPosition( tmpPosition )
    #-------------------------------

    def setExternalBias(self,externalBiasClass):
        self.externalBiasClass = externalBiasClass
        self.useExternalBias = True
    #-------------------------------

    def initPotentialGrid(self):
        self.gridBins = np.full(self.dimension, self.numberOfGridBins, dtype=np.int64)
        self.gridLimitsMin = self.potentialClass.getGridLimitsMin()
        self.gridLimitsMax = self.potentialClass.getGridLimitsMax()
        if self.dimension == 1:
            self.initPotentialGrid1D()
        elif self.dimension == 2:
            self.initPotentialGrid2D()
        else:
            pass
        #-------------------------------

    def initPotentialGrid1D(self):
        self.xGrid = np.linspace(self.gridLimitsMin[0], self.gridLimitsMax[0], self.gridBins[0])
        self.potGrid = np.zeros(self.gridBins)
        for i in range(self.gridBins[0]):
            self.potGrid[i] = self.potentialClass.getPotentialValue([ self.xGrid[i] ])
    #-------------------------------

    def initPotentialGrid2D(self):
        xTmpGrid = np.linspace(self.gridLimitsMin[0], self.gridLimitsMax[0], self.gridBins[0])
        yTmpGrid = np.linspace(self.gridLimitsMin[1], self.gridLimitsMax[1], self.gridBins[1])
        [ self.xGrid, self.yGrid ] = np.meshgrid(xTmpGrid,yTmpGrid)
        self.potGrid = np.zeros(self.gridBins)
        for i in range(self.gridBins[0]):
            for j in range(self.gridBins[1]):
                self.potGrid[i,j] = self.potentialClass.getPotentialValue([ self.xGrid[i,j], self.yGrid[i,j] ])
    #-------------------------------

    def getPotentialGrid(self,minToZero=False):
        outPotGrid = np.copy(self.potGrid)
        if minToZero:
            outPotGrid = outPotGrid - np.min(outPotGrid)
        if self.dimension == 1:
            return (np.copy(self.xGrid), outPotGrid)
        elif self.dimension == 2:
            return (np.copy(self.xGrid), np.copy(self.yGrid), outPotGrid)
        else:
            pass
    #-------------------------------

    def plotPotentialAndTrajectory(self,FigNumber=1,plotTrajectory=True,FigFilename=None):
        if self.dimension == 1:
            self.plotPotentialAndTrajectory1D(FigNumber,plotTrajectory,FigFilename)
        elif self.dimension == 2:
            self.plotPotentialAndTrajectory2D(FigNumber,plotTrajectory,FigFilename)
        else:
            pass
    #-------------------------------

    def plotPotentialAndTrajectory1D(self,FigNumber,plotTrajectory,FigFilename):
        plt.figure(FigNumber)
        plt.clf()
        yLim = self.potentialClass.getPlotLimits()
        plt.plot(self.xGrid,self.potGrid,lw=3, color='black')
        plt.ylim(yLim)
        if plotTrajectory:
            plt.plot(self.trajectory,self.potentialValues,
                marker=self.particleMarker,
                markersize=self.particleMarkerSize,
                markeredgecolor=self.particleMarkerColor,
                markerfacecolor=self.particleMarkerColor
                )
        plt.xlabel('x')
        plt.ylabel('V(x)')
        if FigFilename != None:
            FigFormat = FigFilename[-3:]
            plt.savefig(FigFilename, transparent=True, format=FigFormat)
        plt.show()
    #-------------------------------

    def plotPotentialAndTrajectory2D(self,FigNumber,plotTrajectory,FigFilename):
        plt.figure(FigNumber)
        plt.clf()
        zLim = self.potentialClass.getPlotLimits()
        plt.contourf(self.xGrid, self.yGrid, self.potGrid,
                     vmin=zLim[0], vmax=zLim[1],
                     levels=np.linspace(zLim[0],zLim[1],1000) )
        plt.colorbar(ticks=np.linspace(zLim[0],zLim[1],5), boundaries=np.linspace(-zLim[0],zLim[1],5))
        if plotTrajectory:
            xTraj = self.getTrajectoryArray()[:,0]
            yTraj = self.getTrajectoryArray()[:,1]
            plt.plot( xTraj, yTraj,
                lw=4,
                marker=self.particleMarker,
                markersize=self.particleMarkerSize,
                markeredgecolor=self.particleMarkerColor,
                markerfacecolor=self.particleMarkerColor
                )
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis([self.gridLimitsMin[0], self.gridLimitsMax[0], self.gridLimitsMin[1], self.gridLimitsMax[1]])
        if FigFilename != None:
            FigFormat = FigFilename[-3:]
            plt.savefig(FigFilename, transparent=True, format=FigFormat)
        plt.show()
    #-------------------------------

    def plotTrajectoryTimeSeries(self,FigNumber=10):
        ylabel = ['x','y','z']
        for i in range(self.dimension):
            plt.figure(FigNumber+i)
            plt.clf()
            traj = self.getTrajectoryArray()[:,i]
            plt.plot(traj)
            ylim_min = self.potentialClass.getGridLimitsMin()[i]
            ylim_max = self.potentialClass.getGridLimitsMax()[i]
            plt.ylim(ylim_min,ylim_max)
            plt.xlabel('Time')
            plt.ylabel(ylabel[i]+' coordinate')
        plt.show()
    #-------------------------------

    def getPosition(self):
        return np.copy(self.position)
    #-------------------------------

    def setPosition(self,value):
        assert len(value) == self.dimension
        self.position = np.array(value)
    #-------------------------------

    def getTemperature(self):
        return self.T
    #-------------------------------

    def setTemperature(self,Temperature):
        self.T = float(Temperature)
        self.kBT = self.kB * self.T
        self.beta = 1.0/self.kBT
    #-------------------------------

    def getBeta(self):
        return self.beta
    #-------------------------------

    def getAcceptenceProbability(self,deltaE):
        prob = np.exp(-self.beta*deltaE)
        if prob >= 1.0:
            prob = 1.0
        return prob
    #-------------------------------

    def getStepSize(self):
        return self.MCstepSize
    #-------------------------------

    def setStepSize(self,MCstepSize):
        self.MCstepSize = MCstepSize
    #-------------------------------

    def randomDisplacement(self,oldPosition):
        newPosition = np.zeros(self.dimension)
        displacementRMSD = 0.0
        for i in range(self.dimension):
            delta = self.MCstepSize * np.random.uniform(-1.0,1.0)
            displacementRMSD += delta**2
            newPosition[i] = oldPosition[i] + delta
        displacementRMSD = np.sqrt(displacementRMSD/self.dimension)
        return (newPosition, displacementRMSD)
    #-------------------------------

    def doOneMonteCarloMove(self):
        self.movesTotal += 1
        # get current and trial position
        currPosition = self.getPosition()
        (trialPosition, displacementRMSD) = self.randomDisplacement(currPosition)
        # get values of potential and external bias
        currPot = self.potentialClass.getPotentialValue(currPosition)
        trialPot = self.potentialClass.getPotentialValue(trialPosition)
        currExtBias = 0.0
        trialExtBias = 0.0
        if self.useExternalBias:
            currExtBias = self.externalBiasClass.getBiasValue(currPosition)
            trialExtBias = self.externalBiasClass.getBiasValue(trialPosition)
        # energy difference for the trial move
        deltaPot = (trialPot + trialExtBias) - (currPot + currExtBias)
        # get acceptence probability and update position if accepted
        prob = self.getAcceptenceProbability(deltaPot)
        if prob > np.random.uniform(0.0,1.0):
            self.movesAccepted += 1
            self.setPosition(trialPosition)
            currPot = trialPot
            currExtBias = trialExtBias
        else:
            self.movesRejected += 1
            displacementRMSD = 0.0
        # add to trajectory
        self.trajectory.append( self.getPosition() )
        self.potentialValues.append( currPot )
        self.externalBiasValues.append( currExtBias )
        self.totalDisplacement += displacementRMSD
        # check if the external bias is to be updated
        if self.useExternalBias:
            self.externalBiasClass.updateBias(self.movesTotal,self.getPosition())
    #-------------------------------

    def runMC(self,numberOfMoves,verbose=True):
        for i in range(numberOfMoves):
            self.doOneMonteCarloMove()
            if verbose and (i+1) % (numberOfMoves/10) == 0:
                print "{:8d} of {:8d} steps done".format(i+1,numberOfMoves)
    #-------------------------------

    def getAverageAcceptence(self):
        aveAcceptence = ( np.float(self.movesAccepted) / np.float(self.movesTotal) )
        return aveAcceptence
    #-------------------------------

    def printAverageAcceptence(self):
        AveAcceptence = self.getAverageAcceptence()
        print " Average Acceptence: {:f}".format(AveAcceptence)
    #-------------------------------

    def randomPosition(self):
        randomPosition = np.zeros(self.dimension)
        for i in range(self.dimension):
            min = self.potentialClass.getGridLimitsMin()[i]
            max = self.potentialClass.getGridLimitsMax()[i]
            randomPosition[i] = np.random.uniform(min,max)
        return randomPosition
    #-------------------------------

    def resetPosition(self):
        self.position = np.zeros(self.dimension)
    #-------------------------------

    def resetCountersAndAverages(self):
        self.movesTotal = 0
        self.movesAccepted = 0
        self.movesRejected = 0
        self.totalDisplacement = 0.0
    #-------------------------------

    def resetTrajectory(self):
        self.trajectory = []
        self.potentialValues = []
        self.externalBiasValues = []
    #-------------------------------

    def resetRun(self):
        self.resetCountersAndAverages()
        self.resetTrajectory()
    #-------------------------------

    def getTrajectoryArray(self):
        return np.array(self.trajectory)
    #-------------------------------

    def getTrajectoryHistogramAndFES(self):
        if self.dimension == 1:
            (histoGrid, xEdges) = np.histogram(self.trajectory,
                                    bins=self.gridBins[0],
                                    range=(self.gridLimitsMin[0],self.gridLimitsMax[0]),
                                    density=True)
            old_err_settings = np.seterr(all='ignore')
            fesGrid = - self.kBT * np.log(histoGrid)
            np.seterr(**old_err_settings)
            fesGrid = fesGrid - np.min(fesGrid)
            gridSpacing = xEdges[1]-xEdges[0]
            xGrid = np.linspace(xEdges[0]+0.5*gridSpacing,xEdges[-1]-0.5*gridSpacing,self.gridBins[0])
            return (xGrid, histoGrid, fesGrid)
        elif self.dimension == 2:
            pass
        else:
            pass
    #-------------------------------

    def plotTrajectoryHistogramAndFES(self,FigNumber=10):
        if self.dimension == 1:
            self.plotTrajectoryHistogramAndFES_1D()
        elif self.dimension == 2:
            self.plotTrajectoryHistogramAndFES_2D()
        else:
            pass
    #-------------------------------

    def plotTrajectoryHistogramAndFES_1D(self,FigNumber=10):
        (xGrid1, histoGrid, fesGrid) = self.getTrajectoryHistogramAndFES()
        (xGrid2, potGrid) = self.getPotentialGrid(minToZero=True)
        plt.figure(FigNumber)
        plt.clf()
        plt.plot(xGrid2,potGrid,lw=2,color='black')
        plt.plot(xGrid1,fesGrid,lw=4,color='blue')
        yLim = self.potentialClass.getPlotLimits()
        plt.ylim(yLim)
        plt.show()
    #-------------------------------

    def plotTrajectoryHistogramAndFES_2D(self,FigNumber=10):
        pass
    #-------------------------------

    def getTrajectoryMeanAndStdDev(self):
        traj = self.getTrajectoryArray()
        mean = np.zeros(self.dimension)
        stddev = np.zeros(self.dimension)
        for i in range(self.dimension):
            mean[i] = np.mean(traj[:,i])
            stddev[i] = np.std(traj[:,i])
        return (mean, stddev)
    #-------------------------------

    def printTrajectoryMeanAndStddev(self):
        if self.dimension == 1:
            self.printTrajectoryMeanAndStddev_1D()
        elif self.dimension == 2:
            self.printTrajectoryMeanAndStddev_2D()
        else:
            self.printTrajectoryMeanAndStddev_General()
        #-------------------------------

    def printTrajectoryMeanAndStddev_1D(self):
        (mean, stddev) = self.getTrajectoryMeanAndStdDev()
        print " x position:"
        print "  Mean:   {:15.8f}".format(mean[0])
        print "  StdDev: {:15.8f}".format(stddev[0])
    #-------------------------------

    def printTrajectoryMeanAndStddev_2D(self):
        (mean, stddev) = self.getTrajectoryMeanAndStdDev()
        print " x position:"
        print "  Mean:   {:15.8f}".format(mean[0])
        print "  StdDev: {:15.8f}".format(stddev[0])
        print " y position:"
        print "  Mean:   {:15.8f}".format(mean[1])
        print "  StdDev: {:15.8f}".format(stddev[1])
    #-------------------------------


    def printTrajectoryMeanAndStddev_General(self):
        (mean, stddev) = self.getTrajectoryMeanAndStdDev()
        for i in range(self.dimension):
            print " Dimension {:d}:".format(i+1)
            print "  Mean:   {:15.8f}".format(mean[i])
            print "  StdDev: {:15.8f}".format(stddev[i])
    #-------------------------------

    def getPotentialValuesArray(self):
        return np.array(self.potentialValues)
    #-------------------------------

    def getPotentialValuesMeanAndStddev(self):
        mean = np.mean(self.potentialValues)
        stddev = np.std(self.potentialValues)
        return (mean, stddev)
    #-------------------------------

    def printPotentialValuesMeanAndStddev(self):
        (mean, stddev) = self.getTrajectoryMeanAndStdDev()
        print " Potential Energy:"
        print "  Mean:   {:15.8f}".format(mean)
        print "  StdDev: {:15.8f}".format(stddev)
    #-------------------------------

    def getExternalBiasValuesArray(self):
        return np.array(self.externalBiasValues)
    #-------------------------------

    def getStepNumber(self):
        return self.movesTotal
    #-------------------------------

    def getAverageDisplacement(self):
        value = self.totalDisplacement/float(self.movesTotal)
        return value
    #-------------------------------

    def getPotentialValue(self,position):
        position = np.array(position)
        return self.potentialClass.getPotentialValue(position)
    #-------------------------------

    def getExternalBiasValue(self,position):
        position = np.array(position)
        return self.externalBiasClass.getBiasValue(position)
    #-------------------------------
