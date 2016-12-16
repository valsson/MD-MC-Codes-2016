import numpy as np
from TrajectoryIO import *
import matplotlib.pyplot as plt
from DataTools import writeDataToFile

sigma = 1.0
epsilon = 1.0
m = 1.0
kB = 1.0

r_cutoff = float('inf')
r2_cutoff = r_cutoff**2
ener_shift = 4.0*epsilon*(1.0/r2_cutoff**6-1.0/r2_cutoff**3)

r_walls = 3.0
r2_walls = r_walls**2
k_walls = 100.0

def getLennardJonesPotentialAndForces(p,cell):
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


def getConfiningPotentialAndForces(p):
    forces = np.zeros(p.shape)
    pot = 0.0
    for i in range(nump):
        r2 = np.sum(p[i]**2)
        if r2>r2_walls:
            r = np.sqrt(r2)
            pot += 0.5*k_walls*(r-r_walls)**2
            forces[i] = -k_walls*p[i]*(1.0-r_walls/r)
    return (pot,forces)
#-------------------------------

def getPotentialAndForces(p,cell=None):
    (pot_lj, forces_lj) = getLennardJonesPotentialAndForces(p,cell)
    (pot_w, forces_w) = getConfiningPotentialAndForces(p)
    return (pot_lj+pot_w, forces_lj+forces_w)
#-------------------------------

def getKineticEnergy(v):
    kin = 0.0
    for i in range(v.shape[0]):
        kin += m * np.sum(v[i]**2)
    return kin
#-------------------------------

def randomVeloctites(shape,T,seed=None):
    mean = 0.0
    stddev = np.sqrt( (kB*T)/m )
    if seed is not None: np.random.seed(seed)
    return np.random.normal(mean,stddev,shape)
#-------------------------------

def rescaleVeloctites(v,T,dof):
    setVelocityCenterOfMassToZero(v)
    K = getKineticEnergy(v)
    K_target = 0.5*dof*kB*T
    v *= np.sqrt(K_target/K)
#-------------------------------

def setVelocityCenterOfMassToZero(v):
    v_com = np.average(v,axis=0)
    for i in range(v.shape[0]): v[i] -= v_com
#-------------------------------

def setCenterOfMassToZero(p):
    p_com = np.average(p,axis=0)
    for i in range(p.shape[0]): p[i] -= p_com
#-------------------------------

def getDistancesFromOrigin(p):
    return np.sum(p**2,1)
#-------------------------------





T = 0.8
seed = None
dt = 0.005
num_steps = 100000
#
dim = 2
input_geometry = 'LJ7-2D-initial.xyz'
fn_traj = "traj.xyz"
fn_data = 'results.data'
stride_traj = 10
stride_com = 1
stride_rescale_vel = 1000
#

(postions_in, cell_in) = readPostionsFromFileXYZ(input_geometry,dim)
nump = postions_in.shape[0]
dof = dim*nump - 2.0
p = postions_in
cell = cell_in

v = randomVeloctites(p.shape,T,seed)
rescaleVeloctites(v,T,dof)
header = "time={0}"

setCenterOfMassToZero(p)
setVelocityCenterOfMassToZero(v)

potential_energy = np.zeros(num_steps+1)
kinetic_energy = np.zeros(num_steps+1)
total_energy = np.zeros(num_steps+1)
times = np.arange(num_steps+1)*dt

(pot, F) = getPotentialAndForces(p,cell)

potential_energy[0] = pot
kinetic_energy[0] = getKineticEnergy(v)

writePostionsToFileXYZ(fn_traj,p,['Ar'],header.format(times[0]),cell,False)

for i in range(num_steps):
    p += v*dt+0.5*(F/m)*dt**2
    (pot_new, Fnew) = getPotentialAndForces(p,cell)
    v += (0.5/m)*(Fnew+F)*dt
    if (i+1) % stride_com == 0: setVelocityCenterOfMassToZero(v)
    if (i+1) % stride_rescale_vel == 0: rescaleVeloctites(v,T,dof)
    if (i+1) % stride_traj == 0: writePostionsToFileXYZ(fn_traj,p,['Ar'],header.format(times[i+1]),cell,True)
    kinetic_energy[i+1] = getKineticEnergy(v)
    potential_energy[i+1] = pot_new
    F = Fnew


writePostionsToFileXYZ('postions_final.xyz',p,['Ar'],'',cell,False)
writePostionsToFileXYZ('velocites_final.xyz',v,['Ar'],'',None,False)

total_energy = potential_energy + kinetic_energy
temperature = 2.0*kinetic_energy/(dof*kB)

writeDataToFile(fn_data,
                [times,potential_energy,kinetic_energy,total_energy,temperature],
                ['time','pot','kin','total','T'],
                dataFormat='%15.8f')

plt.figure(3)
plt.plot(times,potential_energy)
plt.figure(4)
plt.plot(times,kinetic_energy)
plt.figure(5)
plt.ylim(np.min(total_energy)*1.1, 0.0)
plt.plot(times,total_energy)
plt.figure(6)
plt.plot(times,temperature)
plt.show()
