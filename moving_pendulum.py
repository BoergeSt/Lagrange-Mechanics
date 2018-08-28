from components import Connector, Point, FixLine, FixCircle, Trolley
from simulation import Simulation

import numpy as np
import sympy as sp

sim = Simulation(movie=False)

#Base = FixLine()
#Base = FixLine(point2 = np.array([1,.01]))
trace = [0.3*sp.sin(2*sim.t),sp.Integer(0)]
Base = FixCircle(1/4,midpoint=trace,moving=True)
T1 = Trolley(Base,np.pi,mass = 10)
C1 = Connector(T1,phi0 = 8/9*np.pi)
P1 = Point(C1)
sim.addObjects([Base,T1,C1,P1])

#C2 = Connector(P1, length=0.5)
#P2 = Point(C2)
#sim.addObjects([C2,P2])



sim.run()


