#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from DataTools import writeDataToFile
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--time-step',dest='time_step',required=False)
parser.add_argument('--output-file',dest='fn_out',required=False)
args = parser.parse_args()

# Parameters of potential
m = 1.0
k = (2.0*np.pi)**2
angular_freq = np.sqrt(k/m)
freq = angular_freq/(2.0*np.pi)
period = 1.0/freq

# MD Parameters
if(args.time_step):
    time_step = np.float64(args.time_step)
else:
    time_step = 0.01*period
if(args.fn_out):
    fn_out = args.fn_out
else:
    fn_out = 'results.data'
showPlots = False
#num_periods = 20
#num_steps = np.int(np.rint( (num_periods*period)/time_step   ))
num_steps = 10000
# initial postion and velocity at t=0
initial_position = 2.0
initial_velocity = 0.0

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


# analytical solution:
phi = np.arctan(-initial_velocity/(initial_position*angular_freq))
amplitude =  initial_position/np.cos(phi)
conserved_energy = getPotentialEnergy(amplitude)

# ----------------------
times = []
positions = []
velocites = []
pot_energies = []
kin_energies = []
tot_energies = []

time = 0.0
curr_position = initial_position
curr_velocity = initial_velocity-0.5*getForce(curr_position)*time_step

for i in range(num_steps):
    if (i+1) % (num_steps/10) == 0:
        print 'MD step {0:6d} of {1:6d}'.format(i+1,num_steps)
    # get force at t
    accleration = getAccleration(curr_position)
    # get new velocites at t+0.5dt
    new_velocity = curr_velocity + accleration*time_step
    # get new positions at t+dt
    new_position = curr_position + new_velocity*time_step
    # get energies at t
    curr_pot_ener = getPotentialEnergy(curr_position)
    curr_kin_ener = getKineticEnergy((new_velocity+curr_velocity)/2.0)
    curr_tot_ener = curr_pot_ener + curr_kin_ener
    #
    times.append( time )
    positions.append( curr_position )
    velocites.append( curr_velocity )
    pot_energies.append( curr_pot_ener )
    kin_energies.append( curr_kin_ener )
    tot_energies.append( curr_tot_ener )
    #
    curr_velocity = new_velocity
    curr_position = new_position
    time += time_step
    #
#----------------------------------------

times = np.array(times)
positions = np.array(positions)
velocites = np.array(velocites)
pot_energies = np.array(pot_energies)
kin_energies = np.array(kin_energies)
tot_energies  = np.array(tot_energies)

positions_analytical = amplitude*np.cos(angular_freq*times+phi)
velocites_analytical = -angular_freq*amplitude*np.sin(angular_freq*times+phi)

writeDataToFile(fn_out,
                [times,positions,velocites,pot_energies,kin_energies,tot_energies,positions_analytical,velocites_analytical],
                ['time','pos','vel','pot_ene','kin_ene','tot_ene','pos_an','vel_an'],
                constantsNames=['time_step','period','amplitude','k','m','phi','conserved_energy'],
                constantsValues=[time_step,period,amplitude,k,m,phi,conserved_energy],
                dataFormat='%15.8f')



if showPlots:
    plt.figure(1)
    plt.plot(times,tot_energies)
    plt.plot(times,pot_energies)
    plt.plot(times,kin_energies)
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
    plt.plot(times,positions_analytical)
    plt.show()

    plt.figure(6)
    plt.plot(times,positions-positions_analytical)
    plt.show()
#
