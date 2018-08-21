# Lagrange-Mechanics
A simple Python program, which allows the automatic symbolic creation of the Lagrange equations for pendulums and similar objects. Furthermore a numerical solver is used in order to approximate the solutions.

## Disclaimer
I did this project out of curiosity and interest. I do not claim that any of the results are correct. Furthermore the numerical methods used in this program are not especially suitable for solving chaotic systems like the double pendulum. 

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


## Used Moduls
* sympy: The main component of the program, since it allows symbolic calculations.
* numpy: Used for some numerical calculations
* scipy: Used for the numerical integration of the resulting ODEs and physical constants
* matplotlib: Used for the display of the solution
* time: Well... in order to measure some times...
* logging: In order to log stuff


## Possible Scenarios and example implementations

* Any kind of kombined mathematical pendulum like a double pendulum with varing string length and point mass. An example is implemented in double_pendulum.py and described in the quick start section.
* A pendulum which is suspended from a trolley which can move freely on a static Line or Circle. An example of this is implemented in moving_pendulum.py
* A pendulum propelled by a FixPoint which is moving on a predefined Path. An example of this is implemented in moving_anchor.py (not fully tested)
* Reasonably well approximations of solid connectors with uniform density. This can be achieved by using a great number of points uniformly distributed on the connector. But since the actual kinetic energy is of second order in case of a single pendulum one can use an order 2 quadrature rule for better approximations. In fact the following 3 point quadrature rule can reproduce a single solid pendulum perfectly.
```
B = FixPoint()
C = Connector(B)
P0 = Point(C4,local = (1+np.sqrt(3/5))/2,mass=5/18)
P1 = Point(C4,local = 1/2,mass=4/9)
P2 = Point(C4,local = (1-np.sqrt(3/5))/2,mass=5/18)
```
* Any combination of the above
