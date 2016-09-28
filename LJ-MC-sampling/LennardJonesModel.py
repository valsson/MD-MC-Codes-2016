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
        self.cell = np.array(cell)
        self.ener_cutoff = 0.0
        self.ener_cutoff = self.getLJpotential(self.r_cutoff)
    #-------------------------------

    def getLJpotential(self,r):
        if r <= self.r_cutoff:
            r = self.sigma/r
            return 4.0*self.epsilon*(r**12-r**6)-self.ener_cutoff
        else:
            return 0.0
    #-------------------------------

    def getEnergy(self,postions):
        N = postions.shape[0]
        energy = 0.0;
        for i in xrange(N):
            for j in xrange(i+1,N):
                r = self.getDistancePBC(postions[i],postions[j])
                if r<self.r_cutoff: energy += self.getLJpotential(r)
        return energy
    #-------------------------------

    def getDistance(self,pos_i,pos_j):
        vec = pos_i - pos_j
        r2 = np.sum(vec**2)
        return np.sqrt(r2)
    #-------------------------------

    def getDistancePBC(self,pos_i,pos_j):
        vec = pos_i - pos_j
        vec = vec - np.floor(vec/self.cell) * self.cell
        r2 = np.sum(vec**2)
        return np.sqrt(r2)
    #-------------------------------


#-------------------------------
