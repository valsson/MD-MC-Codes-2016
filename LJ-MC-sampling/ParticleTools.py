import numpy as np


def randomizePostions(postions,cell,random_seed=None):
    if random_seed is not None: np.random.seed(random_seed)
    for i in range(postions.shape[0]):
        postions[i][0] = np.random.uniform(0,cell[0])
        postions[i][1] = np.random.uniform(0,cell[1])
        postions[i][2] = np.random.uniform(0,cell[2])
#-------------------------------

def randomDisplacement(postions,max_displacement,particle_index,random_seed=None):
    indexes = [particle_index]
    if random_seed is not None: np.random.seed(random_seed)
    for k in indexes:
        postions[k][0] += max_displacement * np.random.uniform(-1,1)
        postions[k][1] += max_displacement * np.random.uniform(-1,1)
        postions[k][2] += max_displacement * np.random.uniform(-1,1)
    #-------------------------------

def writePostionsToFile(filename,postions,particle_names,cell=None,append=True):
    num_particles = postions.shape[0]
    if len(particle_names)==1: particle_names = particle_names*num_particles
    if append:
        f = open(filename,'a')
    else:
        f = open(filename,'w')

    f.write("  {0}\n\n".format(num_particles))
    for i in range(num_particles):
        a = particle_names[i]
        p = postions[i]
        if cell is not None: p = p - np.floor(p/cell) * cell
        out_str = "  {0}   {1:20.9f}  {2:20.9f}  {3:20.9f}\n".format(a,p[0],p[1],p[2])
        f.write(out_str)
    f.close()
    #-------------------------------

def getSquaredDistance(pos_i,pos_j,cell=None):
    vec = pos_i - pos_j
    if cell is not None: vec = vec - np.rint(vec/cell) * cell
    return np.sum(vec**2)
#-------------------------------

def getAllSquaredDistances(postions,r2_cutoff=float('inf'),cell=None):
    distances = []
    num_particles = postions.shape[0]
    for i in range(num_particles):
        for j in range(i+1,num_particles):
            r2 = getSquaredDistance(postions[i],postions[j],cell)
            if r2<r2_cutoff: distances.append(r2)
    return np.array(distances)
#-------------------------------

def getAllSquaredDistancesForOneParticle(postions,index,r2_cutoff=float('inf'),cell=None):
    distances = []
    num_particles = postions.shape[0]
    for i in range(num_particles):
        if i==index: continue
        r2 = getSquaredDistance(postions[index],postions[i],cell)
        if r2<r2_cutoff: distances.append(r2)
    return np.array(distances)
#-------------------------------


#-------------------------------
