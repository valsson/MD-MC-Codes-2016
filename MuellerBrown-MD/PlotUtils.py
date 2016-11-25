import numpy as np
import matplotlib.pyplot as plt

class PlotUtils(object):
    particleMarkerColor   = 'red'
    particleMarker   = 'x'
    particleMarkerSize   = 10.0


    def __init__(self,potential,grid_bins):
        assert len(grid_bins) == 2
        self.plot_min = potential.plot_min
        self.plot_max = potential.plot_max
        xTmpGrid = np.linspace(self.plot_min[0], self.plot_max[0], grid_bins[0])
        yTmpGrid = np.linspace(self.plot_min[1], self.plot_max[1], grid_bins[1])
        [self.xGrid, self.yGrid ] = np.meshgrid(xTmpGrid,yTmpGrid)
        self.potGrid = np.zeros(grid_bins)
        self.xForceGrid = np.zeros(grid_bins)
        self.yForceGrid = np.zeros(grid_bins)
        for i in range(grid_bins[0]):
            for j in range(grid_bins[1]):
                (pot, force) = potential.getPotentialAndForces([ self.xGrid[i,j] , self.yGrid[i,j] ])
                self.potGrid[i,j] = pot
                self.xForceGrid[i,j] = force[0]
                self.yForceGrid[i,j] = force[1]
                #print self.xGrid[i,j] , self.yGrid[i,j], self.potGrid[i,j]
    #-------------------------------

    def plotPotential(self,trajectory=[],FigFilename=None):
        plt.figure(1)
        plt.clf()
        zlim = [self.plot_min[2], self.plot_max[2]]
        plt.contourf(self.xGrid, self.yGrid, self.potGrid,
                vmin=zlim[0], vmax=zlim[1],
                levels=np.linspace(zlim[0],zlim[1],1000) )
        plt.colorbar(ticks=np.linspace(zlim[0],zlim[1],5), boundaries=np.linspace(-zlim[0],zlim[1],5))
        if len(trajectory) > 0:
            xTraj = trajectory[:,0]
            yTraj = trajectory[:,1]
            plt.plot( xTraj, yTraj,
                lw=4,
                marker=self.particleMarker,
                markersize=self.particleMarkerSize,
                markeredgecolor=self.particleMarkerColor,
                markerfacecolor=self.particleMarkerColor)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis([self.plot_min[0], self.plot_max[0], self.plot_min[1], self.plot_max[1]])
        if FigFilename != None:
            FigFormat = FigFilename[-3:]
            plt.savefig(FigFilename, transparent=True, format=FigFormat)
        plt.show()
    #-------------------------------


    def plotForces(self,FigFilename1=None,FigFilename2=None):
        plt.figure(1)
        plt.clf()
        zlim = [-100,100]
        plt.contourf(self.xGrid, self.yGrid, self.xForceGrid,
                vmin=zlim[0], vmax=zlim[1],
                levels=np.linspace(zlim[0],zlim[1],1000) )
        plt.colorbar(ticks=np.linspace(zlim[0],zlim[1],5), boundaries=np.linspace(-zlim[0],zlim[1],5))
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis([self.plot_min[0], self.plot_max[0], self.plot_min[1], self.plot_max[1]])
        if FigFilename1 != None:
            FigFormat = FigFilename1[-3:]
            plt.savefig(FigFilename2, transparent=True, format=FigFormat)
        #
        plt.figure(2)
        plt.clf()
        zlim = [-100,100]
        plt.contourf(self.xGrid, self.yGrid, self.yForceGrid,
                vmin=zlim[0], vmax=zlim[1],
                levels=np.linspace(zlim[0],zlim[1],1000) )
        plt.colorbar(ticks=np.linspace(zlim[0],zlim[1],5), boundaries=np.linspace(-zlim[0],zlim[1],5))
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis([self.plot_min[0], self.plot_max[0], self.plot_min[1], self.plot_max[1]])
        if FigFilename2 != None:
            FigFormat = FigFilename2[-3:]
            plt.savefig(FigFilename2, transparent=True, format=FigFormat)
        plt.show()
    #-------------------------------
