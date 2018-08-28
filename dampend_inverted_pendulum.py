from components import FixPoint, Connector, Point
from simulation import Simulation

import numpy as np
import sympy as sp

sim = Simulation(g=1,subintegrations=20)

r = 0.1
omega = 20

trace = [sp.Integer(0),r*sp.sin(omega*sim.t)]
P1 = FixPoint(moving = True, position = trace)
C1 = Connector(P1,phi0 = np.pi*15/16,dampening = 0.2)
P2 = Point(C1)
sim.addObjects([P1,C1,P2])

sim.run()


