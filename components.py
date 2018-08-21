
import numpy as np
import sympy as sp
import scipy.constants

import matplotlib.pyplot as plt


class Component:
    """Base Class for all Components"""

    def setup(self,i,t):
        return i

    def potential_expr(self):
        return 0

    def kinetic_expr(self):
        return 0
    
    def calculate_ode_functions(self,f,L):
        return f

    def substitude_symbols(self,f):
        return f

    def get_symbol(self):
        return []

    def get_x0(self,x0):
        pass

    def update(self,x):
        pass


    
class FixPoint(Component):
    """A non moving Point object to anchor other objects"""
    def __init__(self,position=np.array([0,0])):
        self.position = position

    def get_position(self):
        return self.position

    def get_position_expr(self):
        return self.position

    def plot(self,ax):
        x,y = self.get_position().T
        ax.plot(x,y,'ok')

class FixLine(Component):
    def __init__(self,point1=np.array([0,0]),point2=np.array([1,0])):
        self.point1 = point1
        self.point2 = point2

    def l2g(self,local):
        return (1-local)*self.point1+local*self.point2

    def l2g_expr(self,local):
        return self.l2g(local)

    def plot(self,ax):
        x,y = np.array([self.l2g(-10),self.l2g(10)]).T
        ax.plot(x,y,'-k')


class FixCircle(Component):
    def __init__(self,radius = 0.5,midpoint=np.array([0,0])):
        self.midpoint = midpoint
        self.radius = radius

    def l2g(self,local):
        return self.midpoint+self.radius*np.array([np.sin(local),-np.cos(local)])

    def l2g_expr(self,local):
        return [self.midpoint[0]+self.radius*sp.sin(local),self.midpoint[1]-self.radius*sp.cos(local)]

    def plot(self,ax):
        circle = plt.Circle(self.midpoint,self.radius,color = "black",fill=False)
        ax.add_artist(circle)



class Point(Component):
    """A mass point stationary on a Connector"""
    def __init__(self,parent,local=1,mass=1):
        self.parent=parent
        self.local = local
        self.mass = mass
        
    def get_position(self):
        return self.parent.l2g(self.local)

    def get_position_expr(self):
        return self.parent.l2g_expr(self.local)

    def plot(self,ax):
        #print("Point at {}".format(self.get_position()))
        x,y = self.get_position().T
        ax.plot(x,y,'or')
    
    def setup(self,i,t):
        self.t = t;
        return i

    def potential_expr(self):
        return self.mass*scipy.constants.g*self.get_position_expr()[1]

    def kinetic_expr(self):
        pos = self.get_position_expr()
        v = [sp.diff(pos[0],self.t),sp.diff(pos[1],self.t)]
        return 0.5*self.mass*(v[0]**2+v[1]**2)

    
class Trolley(Component):
    """A mass point moving on a FixLine"""
    def __init__(self,parent,loc0=0,dloc0=0,mass=1):
        self.parent=parent
        self.mass = mass
        self.loc0 = loc0
        self.dloc0 = dloc0
        self.local = loc0
        
    def setup(self,i,t):
        self.index = i
        self.q=sp.Function('q{}'.format(i))
        self.dq = sp.diff(self.q(t),t)
        self.t = t
        self.Q,self.dQ,self.ddQ = sp.symbols('Q{0} dQ{0} ddQ{0}'.format(self.index))
        return i+1

    def get_position(self):
        return self.parent.l2g(self.local)

    def get_position_expr(self):
        return self.parent.l2g_expr(self.q(self.t))

    def plot(self,ax):
        #print("Point at {}".format(self.get_position()))
        x,y = self.get_position().T
        ax.plot(x,y,'or')

    def potential_expr(self):
        return self.mass*scipy.constants.g*self.get_position_expr()[1]

    def kinetic_expr(self):
        pos = self.get_position_expr()
        v = [sp.diff(pos[0],self.t),sp.diff(pos[1],self.t)]
        return 0.5*self.mass*(v[0]**2+v[1]**2)

    def calculate_ode_functions(self,f,L):
        ODE = sp.diff(sp.diff(L,self.dq),self.t) - sp.diff(L,self.q(self.t)) 
        return f+[ODE]

    def substitude_symbols(self,f):
        for i in range(len(f)):
            f[i] = f[i].subs(sp.diff(self.q(self.t),self.t,2),self.ddQ).subs(sp.diff(self.q(self.t),self.t),self.dQ).subs(self.q(self.t),self.Q)
        return f

    def get_symbol(self):
        return [(self.Q,self.dQ,self.ddQ)]

    def get_x0(self,x0):
        x0.extend([self.loc0,self.dloc0])

    def update(self,x):
        self.local = x[2*self.index]

class Connector(Component):
    def __init__(self,parent,length=1,offset = 0,phi0 = 0, dphi0 = 0):
        self.parent = parent
        self.length = length
        self.offset = offset
        self.phi0 = phi0
        self.dphi0 = dphi0
        self.phi = phi0


    def setup(self,i,t):
        self.index = i
        self.q=sp.Function('q{}'.format(i))
        self.dq = sp.diff(self.q(t),t)
        self.t = t
        self.Q,self.dQ,self.ddQ = sp.symbols('Q{0} dQ{0} ddQ{0}'.format(self.index))
        return i+1

    def l2g(self,local): #local to global
        return self.parent.get_position()+self.length*(self.offset+local)*np.array([np.sin(self.phi),-np.cos(self.phi)])

    def l2g_expr(self,local): #local to global expression
        modifier = self.length*(self.offset+local)
        par_pos = self.parent.get_position_expr()
        return [par_pos[0]+modifier*sp.sin(self.q(self.t)),par_pos[1]-modifier*sp.cos(self.q(self.t))]
    
    def plot(self,ax):
        x,y = np.array([self.l2g(0),self.l2g(1)]).T
        ax.plot(x,y,'-k')

    def calculate_ode_functions(self,f,L):
        ODE = sp.diff(sp.diff(L,self.dq),self.t) - sp.diff(L,self.q(self.t)) 
        return f+[ODE]

    def substitude_symbols(self,f):
        for i in range(len(f)):
            f[i] = f[i].subs(sp.diff(self.q(self.t),self.t,2),self.ddQ).subs(sp.diff(self.q(self.t),self.t),self.dQ).subs(self.q(self.t),self.Q)
        return f

    def get_symbol(self):
        return [(self.Q,self.dQ,self.ddQ)]

    def get_x0(self,x0):
        x0.extend([self.phi0,self.dphi0])

    def update(self,x):
        self.phi = x[2*self.index]
