# Lagrange-Mechanics
A simple Python program, which allows the automatic symbolic creation of the Lagrange equations for pendulums and similar objects. Furthermore a numerical solver is used in order to approximate the solutions.

## Quick Start
We build up the System we want to model using points and connectors. For example in case of the double pendulum we first need a reference point, where we fixate the first pendulum. This is done using the FixPoint Class.
```
Base = FixPoint()
```
After that we add a Connector to the FixPoint and give it an initial displacement or velocity.
```
C1 = Component(Base,phi0=np.pi/4)
```
Then we add a second but movable point to the end of the connector.
```
P1 = Point(C1)
```
We repeat the last two steps once more to add the second pendulum. Finally we need to create a simulation class add all the created components and start the simulation.
```
sim = Simulation()
sim.addObjects([Base, C1, P1, C2, P2])
sim.run()
```



