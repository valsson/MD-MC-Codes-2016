import numpy as np
import copy

from ParticleTools import *

def getLJpotential(r2):
    if r2 <= r2_cutoff:
        inv_r2 = sigma**2/r2
        value = 4.0*epsilon*(inv_r2**6-inv_r2**3)-ener_cutoff
    else:
        value = 0.0
    return value
#-------------------------------

def getEnergyFromDistances(squared_distances):
    energy = 0.0;
    for i in range(squared_distances.shape[0]):
        if squared_distances[i]<r2_cutoff: energy += getLJpotential(squared_distances[i])
    return energy
#-------------------------------

# Parameters:
sigma = 1.0
epsilon = 1.0
InitialSeed = -1
kBoltzman = 1.0
temperature = 0.7
particle_type = "Ar"
beta = 1.0/(kBoltzman*temperature)
r_cutoff = 2.5
r2_cutoff = r_cutoff**2
r2_cutoff = float("inf")
ener_cutoff = 0.0
ener_cutoff = getLJpotential(r2_cutoff)
num_particles = 128
BoxLength = 6.1
cell = np.array([BoxLength, BoxLength, BoxLength])
num_mc_sweeps = 1000
max_displacement = 0.2
nlist_update = 5
#
fn_traj = "traj.xyz"
# ---------------------------------------

r_nlist = r_cutoff + 2*nlist_update*max_displacement
if InitialSeed >= 0: np.random.seed(InitialSeed)


postions = np.zeros([num_particles,3])
randomizePostions(postions,cell)
writePostionsToFile(fn_traj,postions,[particle_type],cell,False)

total_moves = 0
accepted_moves = 0

for i in range(num_mc_sweeps):
    for k in range(num_particles):
        total_moves += 1

        rnd_index = np.random.randint(num_particles)
        trial_postions = np.copy(postions)
        randomDisplacement(trial_postions,max_displacement,rnd_index)

        current_r2 =   getAllSquaredDistancesForOneParticle(postions,rnd_index,r2_cutoff,cell)
        current_energy = getEnergyFromDistances(current_r2)
        trial_r2 =     getAllSquaredDistancesForOneParticle(trial_postions,rnd_index,r2_cutoff,cell)
        trial_energy =   getEnergyFromDistances(trial_r2)
        delta_energy = trial_energy - current_energy
        #
        # current_r2_B =   getAllSquaredDistances(postions,r2_cutoff,cell)
        # current_energy_B = getEnergyFromDistances(current_r2_B)
        # trial_r2_B =     getAllSquaredDistances(trial_postions,r2_cutoff,cell)
        # trial_energy_B = getEnergyFromDistances(trial_r2_B)
        # delta_energy_B = trial_energy_B - current_energy_B

        prob = np.exp(-beta*delta_energy)
        if prob>1.0: prob=1.0
        if prob > np.random.uniform(0.0,1.0):
            accepted_moves += 1
            postions = np.copy(trial_postions)
    #--------------------------
    writePostionsToFile(fn_traj,postions,[particle_type],cell,True)
    if i % (num_mc_sweeps/100) == 0:
        acc_ratio = np.float64(accepted_moves)/np.float64(total_moves)
        print '{0:6d} of {1:6d} MC sweeps done: accepted_ratio = {2:7.5f}'.format(i+1,num_mc_sweeps,acc_ratio)
