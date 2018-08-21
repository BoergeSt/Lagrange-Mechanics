from components import Connector, Point, FixLine, Trolley
from simulation import Simulation

import numpy as np


if __name__=="__main__":
    sim = Simulation(movie=False)
    
    Base = FixLine()
    T1 = Trolley(Base)
    C1 = Connector(T1,phi0=np.pi/10)
    P1 = Point(C1)
    sim.addObjects([Base,T1,C1,P1])

    C2 = Connector(P1, length=0.5)
    P2 = Point(C2)
    sim.addObjects([C2,P2])



    sim.run()


