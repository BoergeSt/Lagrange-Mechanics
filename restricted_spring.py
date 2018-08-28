from components import Spring, Trolley, FixPoint, FixLine, Connector, Point
from simulation import Simulation

import numpy as np


sim = Simulation(movie=False)
    
Base = FixPoint(position = [0,1])
L1 = FixLine()
T1 = Trolley(L1, loc0 = -1,mass = 2)
S1 = Spring(Base, T1, k = 20)
sim.addObjects([Base, L1, T1, S1])

C1 = Connector(T1,length = 0.5)
P1 = Point(C1,mass = .1)
sim.addObjects([C1,P1])

sim.run()