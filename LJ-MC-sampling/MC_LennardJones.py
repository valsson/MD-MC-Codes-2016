import numpy as np

from LennardJonesModel import LennardJonesModel
from ParticleContainer import ParticleContainer



# Parameters:
InitialSeed = -1
kBoltzman = 1.0
temperature = 0.7
beta = 1.0/(kBoltzman*temperature)
r_cutoff = 2.5
num_particles = 128
BoxLength = 10
cell = np.array([BoxLength, BoxLength, BoxLength])
num_mc_sweeps = 1000
max_displacement = 0.5
#
fn_traj = "traj.xyz"
# ---------------------------------------

if InitialSeed >= 0: np.random.seed(InitialSeed)

lj = LennardJonesModel(cell=cell,r_cutoff=r_cutoff)
particles = ParticleContainer(num_particles=num_particles,cell=cell,particle_type="LJ")
particles.randomizePostions()
particles.writePostionsToFile(fn_traj,wrap_pbc=True,append=False)


total_moves = 0
accepted_moves = 0

for i in xrange(num_mc_sweeps):
    for i in xrange(num_particles):
        total_moves += 1

        current_postions = particles.postions
        current_energy = lj.getEnergy(current_postions)

        rnd_index = np.random.randint(num_particles)
        trial_postions = particles.newPostionsFromRandomDisplacement(max_displacement,particle_index=rnd_index)
        trial_energy = lj.getEnergy(trial_postions)

        delta_energy = trial_energy - current_energy
        prob = np.exp(-beta*delta_energy)
        if prob>1.0: prob=1.0
        if prob > np.random.uniform(0.0,1.0):
            accepted_moves += 1
            particles.postion = trial_postions
    #--------------------------
    particles.writePostionsToFile(fn_traj,wrap_pbc=True,append=True)

print np.float64(accepted_moves)/np.float64(total_moves)
