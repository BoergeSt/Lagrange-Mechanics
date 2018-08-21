from components import Connector, Point, FixPoint
from simulation import Simulation

import numpy as np
import sympy as sp


if __name__=="__main__":
    sim = Simulation(movie=False)
    
    trace = [0.3*sp.sin(10*sim.t),sp.Integer(0)]
    Base = FixPoint(moving=True,position=trace)
    C1 = Connector(Base)
    P1 = Point(C1)
    sim.addObjects([Base, C1, P1])

    sim.run()


