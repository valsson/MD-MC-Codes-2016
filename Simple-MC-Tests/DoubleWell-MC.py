#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

num_steps=10000
delta_r = 2.0
T=5.0
beta=1.0/T
initial_position = 0.5

def potential(x):
    return 0.045*x**4-x**2
#---------------------

positions = []
total_moves = 0
accepted_moves = 0

current_position = initial_position

for i in range(num_steps):
    trial_position = current_position + delta_r * np.random.uniform(-1.0,1.0)

    current_potential = potential(current_position)
    trial_potential = potential(trial_position)

    acceptance_prob = np.exp(-beta*( trial_potential - current_potential) )
    if acceptance_prob > 1.0: acceptance_prob = 1.0
    rnd = np.random.uniform(0.0,1.0)
    if acceptance_prob > rnd:
        accepted_moves += 1
        current_position = trial_position
    total_moves += 1
    positions.append(current_position)

aver_acceptance = np.float64(accepted_moves)/np.float64(total_moves)
print 'average acceptance prob = {0}'.format(aver_acceptance)

positions = np.array(positions)

plt.figure(1)
plt.plot(positions)
plt.show()
