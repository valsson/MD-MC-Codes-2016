import numpy as np
import matplotlib.pyplot as plt


# Parameters of potential
m = 1.0
k = 10.0
angular_freq = np.sqrt(k/m)
freq = angular_freq/(2.0*np.pi)
period = 1.0/freq


print angular_freq,period

def getPotentialEnergy(x):
    potential_ener = 0.5*k*x**2
    return potential_ener
#-------------------------------

def getForce(x):
    force = -k*x
    return force
#-------------------------------

def getAccleration(x):
    return getForce(x)/m
#-------------------------------

def getPotentialAndForce(x):
    return ( getPotentialEnergy(x), getForce(x) )
#-------------------------------

def getKineticEnergy(v):
    kinetic_ener = 0.5*m*v**2
    return kinetic_ener
#-------------------------------

def getTotalEnergy(x,v):
    return getPotentialEnergy(x)+getKineticEnergy(v)
#-------------------------------



#
num_periods = 5
time_step = 0.5
# num_steps = 1000
num_steps = np.int(np.rint( (num_periods*period)/time_step   ))
initial_position = 1.0
#initial_prev_position = 2.0
#initial_velocity = (initial_position - initial_prev_position) / time_step
initial_velocity = -100.0
print initial_velocity
# ----------------------

times = []
positions = []
velocites = []
pot_energies = []
kin_energies = []
tot_energies = []


time = 0.0
curr_position = initial_position
# prev_position = initial_prev_position
prev_position = curr_position-initial_velocity*time_step + 0.5*getForce(curr_position)*time_step**2
curr_velocity = initial_velocity

curr_pot_ener = getPotentialEnergy(curr_position)
curr_kin_ener = getKineticEnergy(curr_velocity)
curr_tot_ener = curr_pot_ener + curr_kin_ener

times.append( time )
positions.append( curr_position )
velocites.append( curr_velocity )
pot_energies.append( curr_pot_ener )
kin_energies.append( curr_kin_ener )
tot_energies.append( curr_tot_ener )

for i in range(num_steps):
    if (i+1) % (num_steps/10) == 0:
        print 'MD step {0:6d} of {1:6d}'.format(i+1,num_steps)
    time += time_step
    accleration = getAccleration(curr_position)
    new_position = 2.0*curr_position - prev_position + accleration*time_step**2
    curr_velocity = (new_position - prev_position) / (2.0*time_step)
    prev_position = curr_position
    curr_position = new_position
    #
    curr_pot_ener = getPotentialEnergy(curr_position)
    curr_kin_ener = getKineticEnergy(curr_velocity)
    curr_tot_ener = curr_pot_ener + curr_kin_ener
    print curr_velocity, curr_kin_ener
    #
    times.append( time )
    positions.append( curr_position )
    velocites.append( curr_velocity )
    pot_energies.append( curr_pot_ener )
    kin_energies.append( curr_kin_ener )
    tot_energies.append( curr_tot_ener )
#----------------------------------------

times = np.array(times)
positions = np.array(positions)
velocites = np.array(velocites)
pot_energies = np.array(pot_energies)
kin_energies = np.array(kin_energies)
tot_energies  = np.array(tot_energies)

plt.figure(1)
plt.plot(times,tot_energies)
plt.show()

plt.figure(2)
plt.plot(times,pot_energies)
plt.show()

plt.figure(3)
plt.plot(times,kin_energies)
plt.show()

plt.figure(4)
plt.plot(times,velocites)
plt.show()

plt.figure(5)
plt.plot(times,positions)
plt.show()