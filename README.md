# Lagrange-Mechanics
A simple Python3 program, which allows the automatic symbolic creation of the Lagrange equations for systems containing pendulums, springs and similar objects. It features a very simple way to set up even complex systems. Afterwards the system can be viewed in (close to) realtime, depending on the complexity of the system or recorded and saved as an mp4.

#### Table of Contents
* [Disclamer](#disclaimer)
* [Quick Start](#quick-start)
* [Generation Rules](#generation-rules)
* [Example Systems](#example-systems)
    * [A Single Pendulum](#a-single-pendulum)
    * [The Double Pendulum](#the-double-pendulum)
    * [Pendulum with Moving Base](#pendulum-with-moving-base)
    * [Pendulum Freely Moving on a Circle](#pendulum-freely-moving-on-a-circle)
    * [Dampened Inverted Pendulum](#dampened-inverted-pendulum)
    * [Freely Swinging Spring](#freely-swinging-spring)
    * [Spring Driven Trolly](#spring-driven-trolly)
* [Stuff That May or May Not Be Implemented in the Future](#stuff-that-may-or-may-not-be-implemented-in-the-future)
* [Used Moduls](#used-moduls)


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

## Generation Rules

There are by now a fair number of components from which one can build up a system. Basically they fall into two groups of components (the "FixGroup" and the "PhysicsGroup") which in turn can each be divided into two groups (the "0d" and "1d" objects):

* FixGroup:
  * 0d:
    * FixPoint: Is a point which can be stationary or moving along a predetermined path.
  * 1d:
    * FixLine: An "infinite" Line similar to the FixPoint
    * FixCircle: Like the FixLine but circular. It can also move along predetermined paths.
* PhysicsGroup:
  * 0d:
    * Point: A point which carries a mass. It keeps stationary on the local coordinates of the 1d object it is attached to.
    * Trolley: Like the Point it caries a mass. Additionally it can move freely on the 1d object it is attached to.
  * 1d:
    * Connector: A "finite" line of constant length, which can rotate freely around the point it is connected to.
    * Spring: Like a connector but with variable length, which in turn will create a potential following Hooke's law. It can either have 2 Degrees of freedom (it's rotation and length) if it has one "parent" object or 0 DOFS it is has two parents.

There are some rules for the assembly of a system:
* Components from the FixGroup can not be attached to another object yet other objects can be attached to them.
* Components from the PhysicsGroup must be attached to exactly one other already existing component. Exception: The spring *can* be attached to two different components.
* 1d components always attach to 0d components and vice versa.

Warning:
* Even though you can technically connect a Trolley to a Connector or even Spring, there is currently nothing stoping the trolley from moving along the infinite line spaned by it.

## Example Systems

### A Single Pendulum

The Basis for any System is at least one element of the FixGroup. In this case we only need a point from which the pendulum can be suspended. Then we attache the "String" to this Point which has to be a 1D object from the physics group according to the previously mentioned construction rules. Therefore we can use either the Connector or the spring. Since we do not want the "String" to change length, we choose the Connector. At the end of the Connector we need to put a weight. According to the rules we need to take a 0d physics object. Since we want to stick the mass to the end of the connector and it is supposed to stay there (and not move along the connector), we take the Point. Finally we have to think about the initial conditions under which we want to calculate the system. In this case we want the pendulum to have an initial displacement of pi/4 but no initial velocity.

```python
from components import Connector, Point, FixPoint
from simulation import Simulation

import numpy as np

sim = Simulation()

Base = FixPoint()
C1 = Connector(Base, phi0 = np.pi/4)
P1 = Point(C1)
sim.addObjects([Base, C1, P1])

sim.run()
```

### The Double Pendulum

There is really nothing more to single pendulum except that you can of course change the mass and the length of the components.
```python
from components import Connector, Point, FixPoint
from simulation import Simulation

import numpy as np

sim = Simulation()

Base = FixPoint()
C1 = Connector(Base, phi0 = np.pi/4)
P1 = Point(C1)
sim.addObjects([Base, C1, P1])

C2 = Connector(Base, length = 0.5)
P1 = Point(C2, mass = 2)
sim.addObjects([C2, P2])

sim.run()
```

### Pendulum with Moving Base

This time instead of using an initial displacement or velocity we stimulate the system by moving the base of the pendulum sin-like to the left and right. This is done by setting the moving parameter to True, and using a sympy term with sim.t as variable instead of number for the components of the position.

```python
from components import Connector, Point, FixPoint
from simulation import Simulation

import numpy as np
import sympy as sp

sim = Simulation()

trace = [0.1*sp.sin(3*sim.t),sp.Integer(0)]
Base = FixPoint(moving=True, position=trace)

C1 = Connector(Base)
P1 = Point(C1)

sim.addObjects([Base,C1,P1])
sim.run()

```

### Pendulum Freely Moving on a Circle

This time we don't want the base to be moving in a predetermined way, but to be freely movable on a circle. This can be done using the Trolley in combination with the FixCircl class.

```python
from components import Connector, Point, FixCircle, Trolley
from simulation import Simulation

import numpy as np
import sympy as sp

sim = Simulation()

Base = FixCircle(1/4)
T1 = Trolley(Base)
C1 = Connector(T1,phi0 = np.pi/4)
P1 = Point(C1)
sim.addObjects([Base,T1,C1,P1])

sim.run()

```

### Dampened Inverted Pendulum

Normally a single pendulum has only one stable equilibrium point (which is at the bottom). But if one "shakes" the base of the pendulum fast enough up and down, the upper equilibrium point becomes stable as well. In case of a dampening effect which is generated by using a kind of temporal covariant derivative in the Lagrange equation) we even get an asymptotic equilibrium point. Since the necessary frequency decreases with the local gravity constant, we set it to one for this excample.

```python
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
```

### Freely Swinging Spring

The easiest form of a spring is a spring which can freely swing like a connector but can also vary in length and thereby creating an additional potential. The implementation is nearly identical to the single pendulum:

```python
from components import Spring, Point, FixPoint
from simulation import Simulation

import numpy as np

sim = Simulation()
    
Base = FixPoint()
S1 = Spring(Base, x0=0, phi0 = np.pi/8, k = 100)
P1 = Point(S1)
sim.addObjects([Base,S1,P1])

sim.run()
```

### Spring Driven Trolly

If we just put a Trolly on a horizontal Fixline, nothing is going to happen. However if we add a spring between a point above the line and the Trolly, which is not in its equilibrium state, the spring will pull the Trolly towards itself. Therefore the Trolley will begin to oszillate. In this case, the spring itself has no degrees of freedom of itself.

```python
from components import Spring, Trolley, FixPoint, FixLine, Connector, Point
from simulation import Simulation

import numpy as np

sim = Simulation()
    
Base = FixPoint(position = [0,1])
L1 = FixLine()
T1 = Trolley(L1, loc0 = -1,mass = 2)
S1 = Spring(Base, T1, k = 20)
sim.addObjects([Base, L1, T1, S1])

sim.run()
```

## Stuff That May or May Not Be Implemented in the Future
Likely:
* Connectors with predetermined variing length 
* Using Circles at the end of a connector to which trolleys or points can be connected
* A free particel

Possibly (meaning if I have enough time and motivation... and if I can figure out how to do this):
* Inequality Restraints. For example a trolley that can only move on a finite line which has kind of a stopper at the end. (I don't think anything like a bounce back can be implemented using the lagrange formalism. Yet if you know a way to do that, please let me know.)
* If I become really bored I might add a GUI such that you can assamble your system graphically and change the parameters on the fly.

Maybe (mainly inner workings you probably won't notice):
* Possibility of different general Potential functions.
* Error handling
* Performance improvements such that even complex systems can be viewed in (nearly) realtime


## Used Moduls
* sympy: The main component of the program, since it allows symbolic calculations.
* numpy: Used for some numerical calculations
* scipy: Used for the numerical integration of the resulting ODEs and physical constants
* matplotlib: Used for the display of the solution
* time: Well... in order to measure some times...
* logging: In order to log stuff

Additionally you will need FFmpeg if you want to save the animation as a movie.
