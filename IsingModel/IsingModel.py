import numpy as np
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class IsingModel(object):

    """Ising model"""

    pi = np.pi
    kB = 1.0


    SpinUpColor   = 'blue'
    SpinDownColor = 'red'

    SpinUpMarker   = '^'
    SpinDownMarker = 'v'
    SpinUpMarkerSize   = 100.0
    SpinDownMarkerSize = 100.0
    SpinUpMarkerEdgeColor     = 'none'
    SpinDownMarkerEdgeColor   = 'none'
    SpinUpMarkerFaceColor     = SpinUpColor
    SpinDownMarkerFaceColor   = SpinDownColor
    Scatter3D_Depthshade = True
    Scatter2D_MarkerAlpha = 0.8
    Scatter3D_MarkerAlpha = 0.8

    SpinFigCMAP = ListedColormap([SpinDownColor, SpinUpColor])
    SpinFigLevels = [-1,0,1]



    def __init__(self,
                Dimension = 2,
                L = 21,
                Jcoupling=1.0,
                Hfield = 0.0,
                Temperature = 1.0,
                FigNum = 1,
                InitialState = 'random',
                InitialSeed = -1):
        self.Dim = Dimension
        self.NumNeighb = 2*self.Dim
        self.L = L
        self.NumSpins = L**self.Dim
        self.Jc = float(Jcoupling)
        self.Hfield = float(Hfield)
        self.T = float(Temperature)
        self.beta = 1.0/(self.kB*self.T)
        self.SpinFigNumber = FigNum
        self.shape = []
        for i in range(self.Dim):
            self.shape.append(L)
        self.setupIndices()
        self.setupNearestNeighbourList()
        self.Spins = np.ones(self.NumSpins, dtype=np.int)
        if InitialSeed >= 0: np.random.seed(InitialSeed)

        if InitialState.lower() == 'random':
            self.randomizeSpins()
        elif InitialState.lower() == 'up':
            self.setAllSpinsUp()
        elif InitialState.lower() == 'down':
            self.setAllSpinsDown()
        else:
            pass

        self.MCstep = 0
        self.averagesCounter = 1
        self.initSpinFig()
    #-------------------------------

    def setupIndices(self):
        self.indices2index = np.arange(self.NumSpins, dtype=int).reshape(self.shape)
        self.index2indices = []
        for i in range(self.NumSpins):
            ind = np.array( np.where(self.indices2index == i) ).flatten()
            self.index2indices.append(ind)
        self.index2indices = np.array(self.index2indices)
    #-------------------------------

    def getIndex(self,indices):
        return self.indices2index[indices]
    #-------------------------------

    def getIndicesTuple(self,index):
        return tuple(self.index2indices[index,:])
    #-------------------------------

    def getIndicesArray(self,index):
        return np.copy(self.index2indices[index,:])
    #-------------------------------

    def getNearestNeighbours(self,index):
        NN_indexes = []
        for i in range(self.Dim):
            ind_p = self.getIndicesArray(index)
            ind_m = self.getIndicesArray(index)
            ind_p[i] += 1
            ind_m[i] -= 1
            # Check boundary conditions
            if ind_p[i] == self.L  : ind_p[i] = 0
            if ind_m[i] == -1 : ind_m[i] = self.L-1
            #
            NN_indexes.append( self.getIndex( tuple(ind_p) ) )
            NN_indexes.append( self.getIndex( tuple(ind_m) ) )
        NN_indexes = np.array(NN_indexes)
        NN_indexes.sort()
        return np.copy(NN_indexes)
    #-------------------------------

    def setupNearestNeighbourList(self):
        self.NNlist = []
        for i in range(self.NumSpins):
            self.NNlist.append( self.getNearestNeighbours(i) )
        self.NNlist = np.array(self.NNlist)
    #-------------------------------

    def initSpinFig(self):
        if self.Dim == 1:
            self.initSpinFig1D_Scatter()
        elif self.Dim == 2:
            # self.initSpinFig2D()
            self.initSpinFig2D_Scatter()
        elif self.Dim == 3:
            self.initSpinFig3D_Scatter()
    #-------------------------------

    def initSpinFig1D_Scatter(self):
        self.x = self.index2indices[:,0]
        self.y = np.zeros(self.NumSpins, dtype=np.int)
        indexesUp =  np.where(self.Spins == 1)
        indexesDown = np.where(self.Spins == -1)
        # SpinsUp   = self.Spins[indexesUp]
        # SpinsDown = self.Spins[indexesDown]
        x_up   = self.x[indexesUp]
        x_down = self.x[indexesDown]
        y_up   = self.y[indexesUp]
        y_down = self.y[indexesDown]
        plt.figure(self.SpinFigNumber)
        plt.scatter(x_up, y_up,
                    c=self.SpinUpColor,
                    marker=self.SpinUpMarker,
                    s=self.SpinUpMarkerSize,
                    edgecolors=self.SpinUpMarkerEdgeColor,
                    facecolors=self.SpinUpMarkerFaceColor)
        plt.scatter(x_down, y_down,
                    c=self.SpinDownColor,
                    marker=self.SpinDownMarker,
                    s=self.SpinDownMarkerSize,
                    edgecolors=self.SpinDownMarkerEdgeColor,
                    facecolors=self.SpinDownMarkerFaceColor)
        # plt.scatter(self.x, self.y, c=self.Spins, cmap=self.SpinFigCMAP, vmin=-1, vmax=1)
        # plt.colorbar(ticks=self.SpinFigLevels, boundaries=self.SpinFigLevels)
        # plt.xlabel('x')
        # plt.ylabel('y')
        plt.axis([0-0.5, self.L-1+0.5, 0-0.1, 0+0.1])
        plt.ion()
        plt.show()
    #-------------------------------


    def initSpinFig2D(self):
        [self.X,self.Y] = np.meshgrid( np.arange(self.L, dtype=np.int), np.arange(self.L, dtype=np.int) )
        SpinsTmp = self.Spins.reshape(self.shape)
        plt.figure(self.SpinFigNumber)
        plt.contourf(self.X, self.Y, SpinsTmp,
                     cmap=self.SpinFigCMAP,
                     vmin=-1, vmax=1,
                     levels=self.SpinFigLevels)
        plt.colorbar(ticks=self.SpinFigLevels, boundaries=self.SpinFigLevels)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis([0, self.L-1, 0,self.L-1])
        plt.ion()
        plt.show()
    #-------------------------------

    def initSpinFig2D_Scatter(self):
        self.x = self.index2indices[:,0]
        self.y = self.index2indices[:,1]
        indexesUp = np.where(self.Spins == 1)
        indexesDown = np.where(self.Spins == -1)
        # SpinsUp   = self.Spins[indexesUp]
        # SpinsDown = self.Spins[indexesDown]
        x_up   = self.x[indexesUp]
        x_down = self.x[indexesDown]
        y_up   = self.y[indexesUp]
        y_down = self.y[indexesDown]
        plt.figure(self.SpinFigNumber)
        plt.scatter(x_up, y_up,
                    c=self.SpinUpColor,
                    marker=self.SpinUpMarker,
                    s=self.SpinUpMarkerSize,
                    edgecolors=self.SpinUpMarkerEdgeColor,
                    facecolors=self.SpinUpMarkerFaceColor,
                    alpha=self.Scatter2D_MarkerAlpha)
        plt.scatter(x_down, y_down,
                    c=self.SpinDownColor,
                    marker=self.SpinDownMarker,
                    s=self.SpinDownMarkerSize,
                    edgecolors=self.SpinDownMarkerEdgeColor,
                    facecolors=self.SpinDownMarkerFaceColor,
                    alpha=self.Scatter2D_MarkerAlpha)
        # plt.scatter(self.x, self.y, c=self.Spins, cmap=self.SpinFigCMAP, vmin=-1, vmax=1)
        # plt.colorbar(ticks=self.SpinFigLevels, boundaries=self.SpinFigLevels)
        # plt.xlabel('x')
        # plt.ylabel('y')
        plt.axis([0-0.5, self.L-1+0.5, 0-0.5,self.L-1+0.5])
        plt.ion()
        plt.show()
    #-------------------------------

    def initSpinFig3D_Scatter(self):
        self.x = self.index2indices[:,0]
        self.y = self.index2indices[:,1]
        self.z = self.index2indices[:,2]
        indexesUp = np.where(self.Spins == 1)
        indexesDown = np.where(self.Spins == -1)
        # SpinsUp   = self.Spins[indexesUp]
        # SpinsDown = self.Spins[indexesDown]
        x_up   = self.x[indexesUp]
        x_down = self.x[indexesDown]
        y_up   = self.y[indexesUp]
        y_down = self.y[indexesDown]
        z_up   = self.z[indexesUp]
        z_down = self.z[indexesDown]
        fig = plt.figure(self.SpinFigNumber)
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x_up, y_up, z_up,
                    c=self.SpinUpColor,
                    marker=self.SpinUpMarker,
                    s=self.SpinUpMarkerSize,
                    edgecolors=self.SpinUpMarkerEdgeColor,
                    facecolors=self.SpinUpMarkerFaceColor,
                    alpha=self.Scatter3D_MarkerAlpha,
                    depthshade= self.Scatter3D_Depthshade)
        ax.scatter(x_down, y_down, z_down,
                    c=self.SpinDownColor,
                    marker=self.SpinDownMarker,
                    s=self.SpinDownMarkerSize,
                    edgecolors=self.SpinDownMarkerEdgeColor,
                    facecolors=self.SpinDownMarkerFaceColor,
                    alpha=self.Scatter3D_MarkerAlpha,
                    depthshade= self.Scatter3D_Depthshade)
        # ax.scatter(self.x, self.y, self.z, c=self.Spins, cmap=self.SpinFigCMAP, vmin=-1, vmax=1, s=100.0)
        # plt.axis([0-0.5, self.L-1+0.5, 0-0.5,self.L-1+0.5])
        plt.ion()
        plt.show()
    #-------------------------------

    def updateSpinFig(self):
        if self.Dim == 1:
            self.updateSpinFig1D_Scatter()
        elif self.Dim == 2:
            # self.updateSpinFig2D()
            self.updateSpinFig2D_Scatter()
        elif self.Dim == 3:
            self.updateSpinFig3D_Scatter()
    #-------------------------------

    def updateSpinFig1D_Scatter(self):
        indexesUp = np.where(self.Spins == 1)
        indexesDown = np.where(self.Spins == -1)
        # SpinsUp   = self.Spins[indexesUp]
        # SpinsDown = self.Spins[indexesDown]
        x_up   = self.x[indexesUp]
        x_down = self.x[indexesDown]
        y_up   = self.y[indexesUp]
        y_down = self.y[indexesDown]
        plt.figure(self.SpinFigNumber)
        plt.clf()
        plt.scatter(x_up, y_up,
                    c=self.SpinUpColor,
                    marker=self.SpinUpMarker,
                    s=self.SpinUpMarkerSize,
                    edgecolors=self.SpinUpMarkerEdgeColor,
                    facecolors=self.SpinUpMarkerFaceColor)
        plt.scatter(x_down, y_down,
                    c=self.SpinDownColor,
                    marker=self.SpinDownMarker,
                    s=self.SpinDownMarkerSize,
                    edgecolors=self.SpinDownMarkerEdgeColor,
                    facecolors=self.SpinDownMarkerFaceColor)
        # plt.scatter(self.x, self.y, c=self.Spins, cmap=self.SpinFigCMAP, vmin=-1, vmax=1)
        # plt.colorbar(ticks=self.SpinFigLevels, boundaries=self.SpinFigLevels)
        # plt.xlabel('x')
        # plt.ylabel('y')
        plt.axis([0-0.5, self.L-1+0.5, 0-0.1, 0+0.1])
        plt.draw()
    #-------------------------------


    def updateSpinFig2D(self):
        SpinsTmp = self.Spins.reshape(self.shape)
        plt.figure(self.SpinFigNumber)
        plt.clf()
        plt.contourf(self.X, self.Y, SpinsTmp,
                     cmap=self.SpinFigCMAP,
                     vmin=-1, vmax=1,
                     levels=self.SpinFigLevels)
        plt.colorbar(ticks=self.SpinFigLevels, boundaries=self.SpinFigLevels)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis([0, self.L-1, 0,self.L-1])
        plt.draw()
    #-------------------------------

    def updateSpinFig2D_Scatter(self):
        indexesUp = np.where(self.Spins == 1)
        indexesDown = np.where(self.Spins == -1)
        # SpinsUp   = self.Spins[indexesUp]
        # SpinsDown = self.Spins[indexesDown]
        x_up   = self.x[indexesUp]
        x_down = self.x[indexesDown]
        y_up   = self.y[indexesUp]
        y_down = self.y[indexesDown]
        plt.figure(self.SpinFigNumber)
        plt.clf()
        plt.scatter(x_up, y_up,
                    c=self.SpinUpColor,
                    marker=self.SpinUpMarker,
                    s=self.SpinUpMarkerSize,
                    edgecolors=self.SpinUpMarkerEdgeColor,
                    facecolors=self.SpinUpMarkerFaceColor,
                    alpha=self.Scatter2D_MarkerAlpha)
        plt.scatter(x_down, y_down,
                    c=self.SpinDownColor,
                    marker=self.SpinDownMarker,
                    s=self.SpinDownMarkerSize,
                    edgecolors=self.SpinDownMarkerEdgeColor,
                    facecolors=self.SpinDownMarkerFaceColor,
                    alpha=self.Scatter2D_MarkerAlpha)
        # plt.scatter(self.x, self.y, c=self.Spins, cmap=self.SpinFigCMAP, vmin=-1, vmax=1)
        # plt.colorbar(ticks=self.SpinFigLevels, boundaries=self.SpinFigLevels)
        # plt.xlabel('x')
        # plt.ylabel('y')
        plt.axis([0-0.5, self.L-1+0.5, 0-0.5,self.L-1+0.5])
        plt.draw()
    #-------------------------------

    def updateSpinFig3D_Scatter(self):
        indexesUp = np.where(self.Spins == 1)
        indexesDown = np.where(self.Spins == -1)
        # SpinsUp   = self.Spins[indexesUp]
        # SpinsDown = self.Spins[indexesDown]
        x_up   = self.x[indexesUp]
        x_down = self.x[indexesDown]
        y_up   = self.y[indexesUp]
        y_down = self.y[indexesDown]
        z_up   = self.z[indexesUp]
        z_down = self.z[indexesDown]
        fig = plt.figure(self.SpinFigNumber)
        plt.clf()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x_up, y_up, z_up,
                    c=self.SpinUpColor,
                    marker=self.SpinUpMarker,
                    s=self.SpinUpMarkerSize,
                    edgecolors=self.SpinUpMarkerEdgeColor,
                    facecolors=self.SpinUpMarkerFaceColor,
                    alpha=self.Scatter3D_MarkerAlpha,
                    depthshade= self.Scatter3D_Depthshade)
        ax.scatter(x_down, y_down, z_down,
                    c=self.SpinDownColor,
                    marker=self.SpinDownMarker,
                    s=self.SpinDownMarkerSize,
                    edgecolors=self.SpinDownMarkerEdgeColor,
                    facecolors=self.SpinDownMarkerFaceColor,
                    alpha=self.Scatter3D_MarkerAlpha,
                    depthshade= self.Scatter3D_Depthshade)
        # ax.scatter(self.x, self.y, self.z, c=self.Spins, cmap=self.SpinFigCMAP, vmin=-1, vmax=1, s=100.0)
        plt.draw()
    #-------------------------------

    def saveCurrentSpinFigure(self,filename, format='pdf',updateFigure=False):
        if(updateFigure): self.updateSpinFig()
        plt.figure(self.SpinFigNumber)
        plt.savefig(filename, transparent=True, format=format)
    #-------------------------------

    def SpinsToColors(self):
        SpinsColors = []
        for i in range(self.Spins.size):
            if self.Spins[i] == 1:
                SpinsColors.append(self.SpinUpColor)
            elif self.Spins[i] == -1:
                SpinsColors.append(self.SpinDownColor)
        return SpinsColors
    #-------------------------------

    def flipSpin(self,index):
        self.Spins[index] = -self.Spins[index]
    #-------------------------------

    def flipRandomSpin(self):
        index = np.random.randint(self.numberOfSpins() )
        self.Spins[index] = -self.Spins[index]
    #-------------------------------

    def getSpin(self,index):
        return np.copy(self.Spins[index])
    #-------------------------------

    def setSpin(self,index,value):
        if (value == 1 or value == -1):
            self.Spins[index] = value
    #-------------------------------

    def numberOfSpins(self):
        return self.NumSpins
    #-------------------------------

    def getTemperature(self):
        return self.T
    #-------------------------------

    def getBeta(self):
        return self.beta
    #-------------------------------

    def getTotalMagnetization(self):
        return np.sum(self.Spins)
    #-------------------------------

    def getMagnetizationPerSpin(self):
        return np.float64(self.getTotalMagnetization())/np.float64(self.NumSpins)
    #-------------------------------

    def getAbsoluteMagnetizationPerSpin(self):
        return np.float64(np.abs(self.getTotalMagnetization()))/np.float64(self.NumSpins)
    #-------------------------------

    def getMeanFieldEnergy(self,index):
        NN = np.copy(self.NNlist[index,:])
        energy = -0.5 * self.Jc * self.Spins[index] * np.sum( self.Spins[NN] )
        energy += - self.Hfield * self.Spins[index]
        return energy
    #-------------------------------

    def getEnergy(self):
        energy = 0.0
        energy2 = 0.0
        for i in range(self.NumSpins):
            NN = np.copy(self.NNlist[i,:])
            energy += self.Spins[i]*np.sum( self.Spins[NN] )
            energy2 += self.getMeanFieldEnergy(i)
        energy = -0.5 * self.Jc * energy
        energy += - self.Hfield * self.getTotalMagnetization()
        return energy
    #-------------------------------

    def getEnergyPerSpin(self):
        return np.float64(self.getEnergy())/np.float64(self.NumSpins)
    #-------------------------------

    def getDeltaEnergyForSpinFlip(self,index):
        NN = np.copy(self.NNlist[index,:])
        deltaE = +2.0 * self.Jc * self.Spins[index] * np.sum( self.Spins[NN] )
        deltaE +=  +2.0 * self.Hfield * self.Spins[index]
        return deltaE
    #-------------------------------

    def getProbabilityForSpinFlip(self,deltaE):
        prob = np.exp(-self.beta*deltaE)
        if prob >= 1.0:
            prob = 1.0
        return prob
    #-------------------------------

    def tryToFlipSpin(self,index):
        deltaE = self.getDeltaEnergyForSpinFlip(index)
        prob = self.getProbabilityForSpinFlip(deltaE)
        if prob > np.random.uniform(0.0,1.0):
            self.flipSpin(index)
    #-------------------------------

    def doOneMonteCarloSweep(self):
        for i in range(self.NumSpins):
            randomIndex = np.random.randint(self.NumSpins)
            self.tryToFlipSpin(randomIndex)
        self.MCstep += 1
    #-------------------------------

    def doMonteCarloSweeps(self,NumSweeps):
        for k in range(NumSweeps):
            self.doOneMonteCarloSweep()
    #-------------------------------

    def getSpinsState(self):
        return np.copy(self.Spins)
    #-------------------------------

    def setSpinsState(self,SpinsState):
        if SpinsState.size == self.NumSpins:
            self.Spins = np.copy(SpinsState)
    #-------------------------------

    def setAllSpinsUp(self):
        self.Spins = np.ones(self.NumSpins, dtype=np.int)
    #-------------------------------

    def setAllSpinsDown(self):
        self.Spins = - np.ones(self.NumSpins, dtype=np.int)
    #-------------------------------

    def clear(self):
        self.setAllSpinsUp()
    #-------------------------------

    def randomizeSpins(self):
        for i in range(self.Spins.size):
            random_number = np.random.uniform(-1,1)
            if random_number >= 0.0:
                self.Spins[i] = 1
            else:
                self.Spins[i] = -1
    #-------------------------------
