import numpy as np

A = np.array(  [ -200.0 , -100.0 , -175.0 ,  15.0 ] )
a = np.array(  [   -1.0 ,   -1.0 ,   -6.5 ,   0.7 ] )
b = np.array(  [    0.0 ,    0.0 ,   11.0 ,   0.6 ] )
c = np.array(  [  -10.0 ,  -10.0 ,   -6.5 ,   0.7 ] )
x0 = np.array( [    1.0 ,    0.0 ,   -0.5 ,  -1.0 ] )
y0 = np.array( [    0.0 ,    0.5 ,    1.5 ,   1.0 ] )

pot_shift = +30.33319242243656
scaling_factor = 0.2
A *= scaling_factor


minima = (np.array( [ -0.558,  1.442 ] ),
          np.array( [  0.623,  0.028 ] ),
          np.array( [ -0.050,  0.467 ] ))

saddlePoints = (np.array( [ -0.822,  0.624 ] ),
                np.array( [ -0.212,  0.293 ] ))

xmin = -1.5
ymin = -0.5
zmin =  0.0
xmax =  1.5
ymax =  2.5
zmax = 300 * scaling_factor

plot_min = [ xmin, ymin, zmin ]
plot_max = [ xmax, ymax, zmax ]

def getPotential(position):
    assert len(position) == 2
    pot = 0.0
    x = position[0]
    y = position[1]
    for i in range(4):
        pot += A[i] * np.exp( a[i]*(x-x0[i])**2 + b[i]*(x-x0[i])*(y-y0[i]) + c[i]*(y-y0[i])**2 )
    pot += pot_shift
    return pot
#-------------------------------


def getPotentialAndForces(position):
    assert len(position) == 2
    pot = 0.0
    force = np.array([ 0.0 , 0.0 ])
    x = position[0]
    y = position[1]
    for i in range(4):
        exp_tmp1 = np.exp( a[i]*(x-x0[i])**2 + b[i]*(x-x0[i])*(y-y0[i]) + c[i]*(y-y0[i])**2 )
        pot += A[i] * exp_tmp1
        force[0] += -A[i] * ( 2.0*a[i]*(x-x0[i])+ b[i]*(y-y0[i]) ) * exp_tmp1
        force[1] += -A[i] * ( b[i]*(x-x0[i])+ 2.0*c[i]*(y-y0[i]) ) * exp_tmp1
    pot += pot_shift
    return (pot,force)
#-------------------------------
