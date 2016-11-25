import numpy as np

pot_shift = 0.0

minima = (np.array( [ -1.174,  1.477 ] ),
          np.array( [ -0.831, -1.366 ] ),
          np.array( [  1.124, -1.486 ] ))

maxima = (np.array( [  0.100,  0.050 ] ),)

saddlePoints = (np.array( [ -1.013, -0.036 ] ),
                np.array( [  0.093,  0.174 ] ),
                np.array( [ -0.208, -1.407 ] ))



def getPotential(position):
    assert len(position) == 2
    x = position[0]
    y = position[1]
    pot = x**4 + y**4 - 2.0*x**2 - 4.0*y**2 + x*y + 0.3*x + 0.1*y
    return pot
#-------------------------------


def getPotentialAndForces(position):
    assert len(position) == 2
    pot = 0.0
    force = np.array([ 0.0 , 0.0 ])
    x = position[0]
    y = position[1]
    print x,y
    pot = x**4 + y**4 - 2.0*x**2 - 4.0*y**2 + x*y + 0.3*x + 0.1*y
    force[0] = -(4.0*x**3 - 4.0*x + y + 0.3)
    force[1] = -(4.0*y**3 - 8.0*y + x + 0.1)
    return (pot,force)
#-------------------------------
