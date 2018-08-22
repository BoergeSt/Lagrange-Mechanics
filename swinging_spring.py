from components import Spring, Point, FixPoint
from simulation import Simulation

import numpy as np


if __name__=="__main__":
    sim = Simulation(movie=False)
    
    Base = FixPoint()
    S1 = Spring(Base, x0=0, phi0 = np.pi/8, k = 100)
    P1 = Point(S1)
    sim.addObjects([Base,S1,P1])

    sim.run()


