import numpy as np

from LennardJonesModel import LennardJonesModel
from ParticleContainer import ParticleContainer



# Parameters:
InitialSeed = -1
kBoltzman = 1.0
temperature = 0.7
beta = 1.0/(kBoltzman*temperature)
r_cutoff = 2.5
num_particles = 64
BoxLength = 5.1
cell = np.array([BoxLength, BoxLength, BoxLength])
num_mc_sweeps = 100
max_displacement = 1.0
nlist_update = 5
#
fn_traj = "traj.xyz"
# ---------------------------------------

r_nlist = r_cutoff + 2*nlist_update*max_displacement
if InitialSeed >= 0: np.random.seed(InitialSeed)

lj = LennardJonesModel(cell=cell,r_cutoff=r_cutoff)
particles = ParticleContainer(num_particles=num_particles,cell=cell,particle_type="LJ",r_nlist=r_nlist)
particles.randomizePostions()
particles.writePostionsToFile(fn_traj,wrap_pbc=True,append=False)


total_moves = 0
accepted_moves = 0

for i in xrange(num_mc_sweeps):
    for k in xrange(num_particles):
        total_moves += 1

        current_postions = particles.postions
        current_energy = lj.getEnergyFromPostions(current_postions)

        rnd_index = np.random.randint(num_particles)
        trial_postions = particles.newPostionsFromRandomDisplacement(max_displacement,particle_index=rnd_index)
        trial_energy = lj.getEnergyFromPostions(trial_postions)

        delta_energy = trial_energy - current_energy
        prob = np.exp(-beta*delta_energy)
        if prob>1.0: prob=1.0
        if prob > np.random.uniform(0.0,1.0):
            accepted_moves += 1
            particles.postion = trial_postions
    #--------------------------
    particles.writePostionsToFile(fn_traj,wrap_pbc=True,append=True)
    if i % (num_mc_sweeps/100) == 0:
        acc_ratio = np.float64(accepted_moves)/np.float64(total_moves)
        print '{0} of {1} MC sweeps done: accepted_ratio = {2}'.format(i,num_mc_sweeps,acc_ratio)
