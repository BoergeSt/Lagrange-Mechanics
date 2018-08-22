# Lagrange-Mechanics
A simple Python3 program, which allows the automatic symbolic creation of the Lagrange equations for systems containing pendulums, springs and similar objects. It features a very simple way to set up even complex systems. Afterwards the system can be viewed in (close to) realtime, depending on the complexity of the system or recorded and saved as an mp4.

## Disclaimer
I did this project out of curiosity and interest. I do not claim that any of the results are correct. Furthermore the numerical methods used in this program are not especially suitable for solving chaotic systems like the double pendulum. 

## Quick Start
We build up the system we want to model using points and connectors. For example in case of the double pendulum we first need a reference point, where we fixate the first pendulum. This is done using the FixPoint Class.
```
Base = FixPoint()
```
After that we add a Connector to the FixPoint and give it an initial displacement or velocity.
```
C1 = Connector(Base,phi0=np.pi/4)
```
Then we add a second but movable point to the end of the Connector.
```
P1 = Point(C1)
```
We repeat the last two steps once more to add the second pendulum.
```
C2 = Connector(P1)
P2 = Point(C2)
```
Finally we need to create a simulation class add all the created components and start the simulation.
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

Additionally you will need FFmpeg if you want to save the animation as a movie.


## Possible Scenarios and Example Implementations

* Any kind of kombined mathematical pendulum like a double pendulum with varing string length and point mass. An example is implemented in double_pendulum.py and described in the quick start section.
* A pendulum which is suspended from a trolley which can move freely on a static Line or Circle. (example: moving_pendulum.py)
* A pendulum propelled by a FixPoint or FixCircle which is moving on a predefined Path. (example: moving_anchor.py)
* Ideal springs as Connectors with variable length and therefore 2 gerneralized coordinates. (example: swinging_spring.py)
* Ideal springs as connectors between two components without any own DOFs. For example a trolley on a line or circle which is connected via a spring to another fixpoint which attracts it. (example: restricted_spring.py)
* Reasonably well approximations of solid connectors with uniform density. This can be achieved by using a great number of points uniformly distributed on the connector. But since the actual kinetic energy is of second order in case of a single pendulum one can use an order 2 quadrature rule for better approximations. In fact the following 3 point quadrature rule can reproduce a single solid pendulum perfectly.
```
B = FixPoint()
C = Connector(B)
P0 = Point(C4,local = (1+np.sqrt(3/5))/2,mass=5/18)
P1 = Point(C4,local = 1/2,mass=4/9)
P2 = Point(C4,local = (1-np.sqrt(3/5))/2,mass=5/18)
```
* Any combination of the above

## Stuff That May or May Not Be Implemented in the Future
Likely:
* Connectors with predetermined variing length
* Using Circles at the end of a connector to which trolleys or points can be connected

Possibly (meaning if I have enough time and motivation... and if I can figure out how to do this):
* Energy loss due to friction (Rayleigh dissipation function)
* Inequality Restraints (if I can figure out a way to do that). For example a trolley that can only move on a finite line which has kind of a stopper at the end. (I don't think anything like a bounce back can be implemented using the lagrange formalism. Yet if you know a way to do that, please let me know.)
* If I become really bored I might add a GUI such that you can assamble your system graphically and change the parameters on the fly.

Maybe (mainly inner workings you probably won't notice):
* Possibility of different general Potential functions.
* Error handling
* Performance improvements such that even complex systems can be viewed in (nearly) realtime
