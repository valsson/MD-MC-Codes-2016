import numpy as np


class LennardJonesModel(object):

    """Lennard Jones Model"""

    kB = 1.0
    epsilon = 1.0
    sigma = 1.0

    def __init__(
                self,
                cell,
                r_cutoff = 2.5):
        self.r_cutoff = r_cutoff
        self.r2_cutoff = self.r_cutoff**2
        self.cell = np.array(cell)
        self.ener_cutoff = 0.0
        self.ener_cutoff = self.getLJpotential(self.r2_cutoff)
    #-------------------------------

    def getLJpotential(self,r2):
        if r2 <= self.r2_cutoff:
            inv_r2 = self.sigma**2/r2
            return 4.0*self.epsilon*(inv_r2**6-inv_r2**3)-self.ener_cutoff
        else:
            return 0.0
    #-------------------------------

    def getEnergyFromPostions(self,postions):
        N = postions.shape[0]
        energy = 0.0;
        for i in xrange(N):
            for j in xrange(i+1,N):
                r2 = self.getSquaredDistancePBC(postions[i],postions[j])
                if r2<self.r2_cutoff: energy += self.getLJpotential(r2)
        return energy
    #-------------------------------

    def getEnergyFromDistances(self,squared_distances):
        energy = 0.0;
        for i in xrange(squared_distances.shape[0]):
            if r2<self.r2_cutoff: energy += self.getLJpotential(squared_distances[i])
        return energy
    #-------------------------------

    def getSquaredDistance(self,pos_i,pos_j):
        vec = pos_i - pos_j
        return np.sum(vec**2)
    #-------------------------------

    def getSquaredDistancePBC(self,pos_i,pos_j):
        vec = pos_i - pos_j
        vec = vec - np.floor(vec/self.cell) * self.cell
        return np.sum(vec**2)
    #-------------------------------


#-------------------------------
