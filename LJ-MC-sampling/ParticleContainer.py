import numpy as np


class ParticleContainer(object):

    """"""

    def __init__(
                self,
                num_particles,
                cell,
                particle_type):
        self.num_particles = num_particles
        self.cell = np.array(cell)
        self.postions = np.zeros([num_particles,3])
        self.velocites = np.zeros([num_particles,3])
        self.particle_names = []
        for i in xrange(self.num_particles): self.particle_names.append(particle_type)

    #-------------------------------

    def randomizePostions(self,random_seed=-1):
        if random_seed >= 0: np.random.seed(random_seed)
        for i in xrange(self.num_particles):
            self.postions[i][0] = np.random.uniform(0,self.cell[0])
            self.postions[i][1] = np.random.uniform(0,self.cell[1])
            self.postions[i][2] = np.random.uniform(0,self.cell[2])
    #-------------------------------

    def newPostionsFromRandomDisplacement(self,max_displacement,particle_index=-1,random_seed=-1):
        if particle_index < 0:
            indexes = range(self.num_particles)
        else:
            indexes = [particle_index]
        if random_seed >= 0: np.random.seed(random_seed)
        new_postions = self.postions
        for k in indexes:
            new_postions[k][0] += max_displacement * np.random.uniform(-1,1)
            new_postions[k][1] += max_displacement * np.random.uniform(-1,1)
            new_postions[k][2] += max_displacement * np.random.uniform(-1,1)
        return new_postions
    #-------------------------------

    def writePostionsToFile(self,filename,wrap_pbc=True,append=True):
        if append:
            f = open(filename,'a')
        else:
            f = open(filename,'w')
        f.write("  {0}\n\n".format(self.num_particles))
        for i in xrange(self.num_particles):
            a = self.particle_names[i]
            p = self.postions[i]
            if wrap_pbc: p = p - np.floor(p/self.cell) * self.cell
            out_str = "  {0}   {1:20.9f}  {2:20.9f}  {3:20.9f}\n".format(a,p[0],p[1],p[2])
            f.write(out_str)
        f.close()
    #-------------------------------

#-------------------------------
