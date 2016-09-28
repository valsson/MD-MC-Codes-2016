import numpy as np

def getLJpotential(r):
    if r <= rcutoff:
        r = 1.0/r
        return 4.0*epsilon*(r**12-r**6)-ecutoff
    else:
        return 0.0
#-------------------------------

def getEnergy(postions,cell):
    N = postions.shape[0]
    energy = 0.0;
    counter = 0
    for i in xrange(N):
        for j in xrange(i+1,N):
            r = getDistancePBC(i,j,cell)
            energy += getLJpotential(r)
    return energy
#-------------------------------


def getDistance(i,j):
    vec = postions[i]-postions[j]
    r2 = np.sum(vec**2)
    return np.sqrt(r2)
#-------------------------------

def getDistancePBC(i,j,cell):
    vec = postions[i]-postions[j]
    vec = vec - np.floor(vec/cell) * cell
    r2 = np.sum(vec**2)
    return np.sqrt(r2)
#-------------------------------

def randomizePostions(postions,cell):
    for i in xrange(postions.shape[0]):
        postions[i][0] = np.random.uniform(0,cell[0])
        postions[i][1] = np.random.uniform(0,cell[1])
        postions[i][2] = np.random.uniform(0,cell[2])
#-------------------------------

def randomDisplacements(postions,max_displacement):
    for i in xrange(postions.shape[0]):
        postions[i][0] += max_displacement * np.random.uniform(-1,1)
        postions[i][1] += max_displacement * np.random.uniform(-1,1)
        postions[i][2] += max_displacement * np.random.uniform(-1,1)
#-------------------------------



def writePostionsToFile(filename,postions,atomnames,append=True):
    N = postions.shape[0]
    if append:
        f = open(filename,'a')
    else:
        f = open(filename,'w')
    f.write("  {0}\n\n".format(N))
    for i in xrange(N):
        a = atomnames[i]
        p = postions[i]
        out_str = "  {0}   {1:20.9f}  {2:20.9f}  {3:20.9f}\n".format(a,p[0],p[1],p[2])
        f.write(out_str)
    f.close()
#-------------------------------


def writePostionsToFilePBC(filename,postions,atomnames,cell,append=True):
    N = postions.shape[0]
    if append:
        f = open(filename,'a')
    else:
        f = open(filename,'w')
    f.write("  {0}\n\n".format(N))
    for i in xrange(N):
        a = atomnames[i]
        p = postions[i]
        # apply PBC
        p = p - np.floor(p/cell) * cell
        out_str = "  {0}   {1:20.9f}  {2:20.9f}  {3:20.9f}\n".format(a,p[0],p[1],p[2])
        f.write(out_str)
    f.close()
#-------------------------------




# Parameters:
pi = np.pi
kB = 1.0
epsilon = 1.0
sigma = 1.0
AtomType = "LJ"
#
InitialSeed = -1
#
T = 0.7
rcutoff = 2.5
ecutoff = 0.0
ecutoff = getLJpotential(rcutoff)
num_particles = 10
BoxLength = 5
cell = [BoxLength, BoxLength, BoxLength]
num_mc_sweeps = 1000
#
fn_traj = "traj.xyz"
# ---------------------------------------

cell = np.array(cell)
if InitialSeed >= 0: np.random.seed(InitialSeed)
postions = np.zeros([num_particles,3])
atomnames = []
for i in xrange(num_particles): atomnames.append(AtomType)
randomizePostions(postions,cell)

# for i in xrange(num_mc_sweeps):
#     for k in xrange(num_particles):
#         index = np.random.randint(num_particles)
#     randomDisplacements(postions,1.0)
#     print getEnergy(postions,cell)
