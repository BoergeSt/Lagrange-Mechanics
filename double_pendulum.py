from components import Connector, Point, FixPoint
from simulation import Simulation

import numpy as np


sim = Simulation(movie=False)
    
Base = FixPoint()
C1 = Connector(Base, phi0 = np.pi/4)
P1 = Point(C1)
sim.addObjects([Base, C1, P1])

C2 = Connector(P1, length = 0.5)
P2 = Point(C2)
sim.addObjects([C2, P2])

#C3 = Connector(P2,1/4)
#P3 = Point(C3)
#sim.addObjects([C3,P3])

#C4 = Connector(Base, phi0 = np.pi/4) # Connector with uniform mass
#PL0 = Point(C4,local = (1+np.sqrt(3/5))/2,mass=5/18)
#PL1 = Point(C4,local = 1/2,mass=4/9)
#PL2 = Point(C4,local = (1-np.sqrt(3/5))/2,mass=5/18)
#sim.addObjects([C4,PL0,PL1,PL2])

sim.run()


