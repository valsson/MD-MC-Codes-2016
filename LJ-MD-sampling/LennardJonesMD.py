import numpy as np
from TrajectoryIO import *
import matplotlib.pyplot as plt

sigma = 1.0
epsilon = 1.0
m = 1.0
kB = 1.0

r_cutoff = float('inf')
r2_cutoff = r_cutoff**2
ener_shift = 4.0*epsilon*(1.0/r2_cutoff**6-1.0/r2_cutoff**3)

def getPotentialAndForces(p,cell):
    nump = p.shape[0]
    forces = np.zeros(p.shape)
    pot = 0.0
    for i in range(nump):
        for j in range(i+1,nump):
            vec = p[i]-p[j]
            if cell is not None: vec = vec - np.rint(vec/cell) * cell
            r2 = np.sum(vec**2)
            if r2 <= r2_cutoff:
                inv_r2 = sigma**2/r2
                pot += 4.0*epsilon*(inv_r2**6-inv_r2**3)-ener_shift
                ff = 24.0*inv_r2**4*(2.0*inv_r2**3-1.0)
                forces[i] += vec*ff
                forces[j] -= vec*ff
    return (pot,forces)
#-------------------------------

def getKineticEnergy(v):
    kin = 0.0
    for i in range(v.shape[0]):
        kin += m * np.sum(v[i]**2)
    return kin
#-------------------------------

def initializeVeloctites(v,T,seed):
    mean = 0.0
    stddev = np.sqrt( (kB*T)/m )
    np.random.seed(seed)
    v = np.random.normal(mean,stddev,v.shape)
#-------------------------------

T = 0.1
seed = 15234
#
dim = 2
input_geo = 'LJ7-2D-initial.xyz'
fn_traj = "traj.xyz"
#
dt = 0.005
num_steps = 10000

(postions_in, cell_in) = readPostionsFromFileXYZ(input_geo,dim)
nump = postions_in.shape[0]
dof = dim*nump
p = postions_in
cell = cell_in
v = np.zeros(p.shape)
initializeVeloctites(v,T,seed)
header_gro = "LJ" + str(nump) + ": time={0}"

potential_energy = np.zeros(num_steps+1)
kinetic_energy = np.zeros(num_steps+1)
total_energy = np.zeros(num_steps+1)
times = np.arange(num_steps+1)*dt

(pot, F) = getPotentialAndForces(p,cell)

potential_energy[0] = pot
kinetic_energy[0] = getKineticEnergy(v)

writePostionsToFileXYZ(fn_traj,p,['Ar'],None,False)
writePostionsToFileXYZ('forces.data',F,['Ar'],None,False)


for i in range(num_steps):
    p += v*dt+0.5*(F/m)*dt**2
    (pot_new, Fnew) = getPotentialAndForces(p,cell)
    v += (0.5/m)*(Fnew+F)*dt
    kinetic_energy[i+1] = getKineticEnergy(v)
    potential_energy[i+1] = pot_new
    F = Fnew
    writePostionsToFileXYZ(fn_traj,p,['Ar'])

total_energy = potential_energy + kinetic_energy
temperature = kinetic_energy/(dof*kB)

plt.figure(3)
plt.plot(times,potential_energy)
plt.figure(4)
plt.plot(times,kinetic_energy)
plt.figure(5)
#plt.ylim(np.min(total_energy)*1.1, np.max(total_energy)*1.1)
plt.plot(times,total_energy)
plt.figure(6)
plt.plot(times,temperature)
plt.show()
