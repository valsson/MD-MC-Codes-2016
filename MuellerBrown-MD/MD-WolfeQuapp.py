from WolfeQuapp import getPotentialAndForces
from PlotUtils import PlotUtils
import numpy as np
import matplotlib.pyplot as plt
import WolfeQuapp as wqpot


m=1.0

def getKineticEnergy(velocity):
    return 0.5*m*(velocity[0]**2+velocity[1]**2)


dt = 0.01
num_steps = 1000

#initial_position = np.array( [  0.0 ,   0.0 ] )
initial_position = wqpot.saddlePoints[1]
initial_velocity = np.array( [  -0.1  ,  0.1 ] )

position = np.zeros([num_steps+1,2])
velocity = np.zeros([num_steps+1,2])
potential_energy = np.zeros(num_steps+1)
kinetic_energy = np.zeros(num_steps+1)
total_energy = np.zeros(num_steps+1)
times = np.arange(num_steps+1)*dt

time = 0.0
position[0,:] = initial_position
velocity[0,:] = initial_velocity
kinetic_energy[0] = getKineticEnergy(initial_velocity)
(pot, force) = getPotentialAndForces(initial_position)
potential_energy[0] = pot

for i in range(0,num_steps):
    # get position at t+dt
    position[i+1] = position[i] + velocity[i]*dt+0.5*(force/m)*dt**2
    # get velocity at t+dt
    (new_pot, new_force) = getPotentialAndForces(position[i+1])
    velocity[i+1] = velocity[i] + (0.5/m)*(new_force+force)*dt
    # add stuff
    kinetic_energy[i+1] = getKineticEnergy(velocity[i+1])
    potential_energy[i+1] = new_pot
    force = new_force

total_energy = potential_energy + kinetic_energy


pu = PlotUtils(wqpot,[200,200])
pu.plotPotential(trajectory=position)

plt.figure(1)
plt.plot(times,position[:,0])
plt.figure(2)
plt.plot(times,position[:,1])
plt.figure(3)
plt.plot(times,potential_energy)
plt.figure(4)
plt.plot(times,kinetic_energy)
plt.figure(5)
plt.ylim(np.min(total_energy)-1.0, np.max(total_energy)+1.0)
plt.plot(times,total_energy)
plt.show()
