
import numpy as np
import sympy as sp
import scipy.constants

import matplotlib.pyplot as plt


class Component:
    """Base Class for all Components"""

    def setup(self, i, t):
        """The setup function sets up the Component for later use in the simulation. Must be called before any other function.
        
        The setup enumerates the degrees of freedom such that no collisions between them occur in symbolic calculations. It is called from the Simulation::setup function.
        """
        
        self.t = t
        return i

    def potential_expr(self,g):
        """Gives back the potential energy expression for the current component."""
        
        return sp.Integer(0)

    def kinetic_expr(self):
        """Gives back the kinetic energy for the current component."""

        return sp.Integer(0)
    
    def calculate_ode_functions(self, f, L):
        """Calculates the second order odes corresponding to the inner degrees of freedom from the Lagrange Function and adds them to the f array."""

        pass

    def substitude_symbols(self, f):
        """Exchanges the implicit time dependency of the DOFs by symbols"""
        pass

    def get_symbol(self):
        """Returns all used symbols corresponding to this Component"""
        return []

    def get_x0(self, x0):
        """fills the array x0 with the initial conditions"""
        pass

    def update(self, x):
        """updates the current values of the degrees of freedom corresponding to this component form the array x"""
        pass

    def init_plot(self,ax):
        self.ax = ax

    def evaluate(self,L):
        return L 


    
class FixPoint(Component):
    """A non moving Point object to anchor other objects"""
    def __init__(self, position=[0,0], moving = False):
        self.position = np.array(position)
        self.moving = moving

    def get_position(self, t=0):
        if self.moving:
            return np.array([self.position[0].subs(self.t,t),self.position[1].subs(self.t,t)])
        return self.position

    def get_position_expr(self):
        return self.position

    def plot(self, t=0):
        x,y = self.get_position(t).T
        self.plt_data.set_data(x,y)
        return self.plt_data

    def init_plot(self,ax):
        self.ax = ax
        self.plt_data, = ax.plot([],[],'ok')

class FixLine(Component):
    """A fixed linear track on which a Trolley can move"""
    def __init__(self, point1=[0,0], point2=[1,0],moving = (False, False)):
        self.point1 = np.array(point1)
        self.point2 = np.array(point2)
        self.moving = moving

    def l2g(self, local, t=0):
        point1 = self.point1
        point2 = self.point2
        if self.moving[0]:
            point1 = [point1[i].subs(self.t,t) for i in [0,1]]
        if self.moving[1]:
            point2 = [point2[i].subs(self.t,t) for i in [0,1]]
        return (1-local)*np.array(point1)+local*np.array(point2)

    def l2g_expr(self, local):
        return [(1-local)*self.point1[0]+local*self.point2[0],(1-local)*self.point1[1]+local*self.point2[1]]

    def plot(self, t=0):
        x,y = np.array([self.l2g(-10,t), self.l2g(10,t)]).T
        self.plt_data.set_data(x,y)
        return self.plt_data

    def init_plot(self,ax):
        self.ax = ax
        self.plt_data, = ax.plot([],[],'-k')

class FixCurve(Component):
    """A fixed track which follows an arbitrary curve on which a Trolley can move"""
    def __init__(self, curve, var, moving = False,plot_interval = (-2,2),plot_points = 10):
        self.curve = curve
        self.var = var
        self.moving = moving
        self.plot_interval = plot_interval
        self.plot_points = plot_points

    def l2g(self, local, t=0):
        pos = [self.curve[0].subs(self.var,local),self.curve[1].subs(self.var,local)]
        if self.moving:
            pos = [pos[0].subs(self.t,t),pos[1].subs(self.t,t)]
        return np.array(pos)

    def l2g_expr(self, local):
        pos = [self.curve[0].subs(self.var,local),self.curve[1].subs(self.var,local)]
        return pos

    def plot(self, t=0):
        if self.first and not self.moving:
            return self.plt_data
        eval = np.linspace(self.plot_interval[0],self.plot_interval[1],self.plot_points)
        x,y = np.array([self.l2g(i,t) for i in eval]).T
        self.plt_data.set_data(x,y)
        self.first = True
        return self.plt_data

    def init_plot(self,ax):
        self.ax = ax
        self.plt_data, = ax.plot([],[],'-k')
        self.first = False

class FixCircle(Component):
    """A fixed circular crack on which a Trolley can move"""
    def __init__(self, radius = 0.5, midpoint=np.array([0,0]), moving = False):
        self.midpoint = midpoint
        self.radius = radius
        self.moving = moving

    def l2g(self, local, t=0):
        if self.moving:
            return np.array([self.midpoint[0].subs(self.t,t)+self.radius*sp.sin(local),self.midpoint[1].subs(self.t,t)-self.radius*sp.cos(local)])
        return self.midpoint+self.radius*np.array([np.sin(local),-np.cos(local)])

    def l2g_expr(self, local):
        return [self.midpoint[0]+self.radius*sp.sin(local), self.midpoint[1]-self.radius*sp.cos(local)]

    def init_plot(self,ax):
        self.ax = ax
        t = 0
        if self.moving:
            self.circle = plt.Circle((self.midpoint[0].subs(self.t,t), self.midpoint[1].subs(self.t,t)), self.radius,color = "black",fill=False)
        else:
            self.circle = plt.Circle(self.midpoint,self.radius,color = "black",fill=False)
        
        self.ax.add_artist(self.circle)

    def plot(self, t=0):
        x,y = self.midpoint
        if self.moving:
            x = self.midpoint[0].subs(self.t,t)
            y = self.midpoint[1].subs(self.t,t)
        self.circle.center = (x,y)
        return self.circle



class Point(Component):
    """A mass point stationary on a Connector"""
    def __init__(self, parent, local=1, mass=1):
        self.parent=parent
        self.local = local
        self.mass = mass
        
    def get_position(self, t=0):
        return self.parent.l2g(self.local,t)

    def get_position_expr(self):
        return self.parent.l2g_expr(self.local)
    
    def setup(self, i, t):
        self.t = t
        return i

    def potential_expr(self,g):
        return self.mass*g*self.get_position_expr()[1]

    def kinetic_expr(self):
        pos = self.get_position_expr()
        v = [sp.diff(pos[0], self.t), sp.diff(pos[1],self.t)]
        return 0.5*self.mass*(v[0]**2+v[1]**2)
   
   
    def plot(self, t=0):
        x,y = self.get_position(t).T
        self.plt_data.set_data(x,y)
        return self.plt_data

    def init_plot(self,ax):
        self.ax = ax
        self.plt_data, = ax.plot([],[],'og')
    
class Trolley(Component):
    """A mass point moving on a FixLine"""
    def __init__(self, parent, loc0=0, dloc0=0, mass=1):
        self.parent=parent
        self.mass = mass
        self.loc0 = loc0
        self.dloc0 = dloc0
        self.local = loc0
        self.dlocal = dloc0
        
    def setup(self, i, t):
        self.index = i
        self.q=sp.Function('q{}'.format(i))
        self.dq = sp.diff(self.q(t),t)
        self.t = t
        self.Q,self.dQ,self.ddQ = sp.symbols('Q{0} dQ{0} ddQ{0}'.format(self.index))
        return i+1

    def get_position(self, t=0):
        return self.parent.l2g(self.local,t)

    def get_position_expr(self):
        return self.parent.l2g_expr(self.q(self.t))


    def potential_expr(self,g):
        return self.mass*g*self.get_position_expr()[1]

    def kinetic_expr(self):
        pos = self.get_position_expr()
        v = [sp.diff(pos[0],self.t),sp.diff(pos[1],self.t)]
        return 0.5*self.mass*(v[0]**2+v[1]**2)

    def calculate_ode_functions(self, f, L):
        ODE = sp.diff(sp.diff(L,self.dq),self.t) - sp.diff(L,self.q(self.t)) 
        f.append(ODE)

    def substitude_symbols(self, f):
        for i in range(len(f)):
            f[i] = f[i].subs(sp.diff(self.q(self.t),self.t,2),self.ddQ).subs(sp.diff(self.q(self.t),self.t),self.dQ).subs(self.q(self.t),self.Q)
        

    def get_symbol(self):
        return [(self.Q, self.dQ, self.ddQ)]

    def get_x0(self, x0):
        x0.extend([self.loc0, self.dloc0])

    def update(self, x):
        self.local = x[2*self.index]
        self.dlocal = x[2*self.index+1]


    def plot(self, t=0):
        x,y = self.get_position(t).T
        self.plt_data.set_data(x,y)
        return self.plt_data

    def init_plot(self,ax):
        self.ax = ax
        self.plt_data, = ax.plot([],[],'or')

    def evaluate(self,L):
        return L.subs(self.dQ,self.dlocal).subs(self.dq,self.dlocal).subs(self.Q,self.local).subs(self.q(self.t),self.local)


class Connector(Component):
    def __init__(self, parent, length=1, offset = 0,phi0 = 0, dphi0 = 0, dampening = 0):
        self.parent = parent
        self.length = length
        self.offset = offset
        self.phi0 = phi0
        self.dphi0 = dphi0
        self.phi = phi0
        self.dphi = dphi0
        self.dampening = dampening


    def setup(self, i, t):
        self.index = i
        self.q=sp.Function('q{}'.format(i))
        self.dq = sp.diff(self.q(t),t)
        self.t = t
        self.Q,self.dQ,self.ddQ = sp.symbols('Q{0} dQ{0} ddQ{0}'.format(self.index))
        return i+1

    def l2g(self, local, t=0): #local to global
        return self.parent.get_position(t)+self.length*(self.offset+local)*np.array([np.sin(self.phi),-np.cos(self.phi)])

    def l2g_expr(self, local): #local to global expression
        modifier = self.length*(self.offset+local)
        par_pos = self.parent.get_position_expr()
        return [par_pos[0]+modifier*sp.sin(self.q(self.t)),par_pos[1]-modifier*sp.cos(self.q(self.t))]
    

    def calculate_ode_functions(self,f,L):
        ODE = sp.diff(sp.diff(L,self.dq),self.t) + self.dampening*sp.diff(L,self.dq) - sp.diff(L,self.q(self.t)) 
        f.append(ODE)

    def substitude_symbols(self, f):
        for i in range(len(f)):
            f[i] = f[i].subs(sp.diff(self.q(self.t),self.t,2),self.ddQ).subs(sp.diff(self.q(self.t),self.t),self.dQ).subs(self.q(self.t),self.Q)
        

    def get_symbol(self):
        return [(self.Q, self.dQ, self.ddQ)]

    def get_x0(self, x0):
        x0.extend([self.phi0, self.dphi0])

    def update(self, x):
        self.phi = x[2*self.index]
        self.dphi = x[2*self.index+1]

    def plot(self, t=0):
        x,y = np.array([self.l2g(0,t),self.l2g(1,t)]).T
        self.plt_data.set_data(x,y)
        return self.plt_data

    def init_plot(self,ax):
        self.ax = ax
        self.plt_data, = ax.plot([],[],'-k')

    def evaluate(self,L):
        return L.subs(self.dQ,self.dphi).subs(self.dq,self.dphi).subs(self.Q,self.phi).subs(self.q(self.t),self.phi)


class Spring(Component):
    def __init__(self, parent, secondary_parent = False, length = 1, k=1, x0 = 0, dx0 = 0, phi0 = 0, dphi0 = 0):
        self.parent = parent
        self.secondary_parent = secondary_parent 
        self.length = length
        self.x0 = x0 
        self.dx0 = dx0
        self.x = x0
        self.phi0 = phi0
        self.dphi0 = dphi0
        self.phi = phi0
        self.dphi = dphi0
        self.k = k

    def setup(self, i, t):
        self.t = t
        if not self.secondary_parent:
            self.x_index = i
            self.q1=sp.Function('q{}'.format(self.x_index))
            self.dq1 = sp.diff(self.q1(t),t)
            self.Q1,self.dQ1,self.ddQ1 = sp.symbols('Q{0} dQ{0} ddQ{0}'.format(self.x_index))
            i += 1;

            self.phi_index = i
            self.q2=sp.Function('q{}'.format(self.phi_index))
            self.dq2 = sp.diff(self.q2(t),t)
            self.Q2,self.dQ2,self.ddQ2 = sp.symbols('Q{0} dQ{0} ddQ{0}'.format(self.phi_index))
            i += 1
        return i

    def l2g(self, local, t=0):
        if not self.secondary_parent:
            return self.parent.get_position(t)+(self.length+self.x)*(local)*np.array([np.sin(self.phi),-np.cos(self.phi)])
        return (1-local)*self.parent.get_position(t)+local*self.secondary_parent.get_position(t)

    def l2g_expr(self, local):
        par_pos = self.parent.get_position_expr()
        if not self.secondary_parent:
            modifier = (self.length+self.q1(self.t))*local
            return [par_pos[0]+modifier*sp.sin(self.q2(self.t)),par_pos[1]-modifier*sp.cos(self.q2(self.t))]
        par_pos2 = self.secondary_parent.get_position_expr()
        return [(1-local)*par_pos[0]+local*par_pos2[0], (1-local)*par_pos[1]+local*par_pos2[1]]
    



    def calculate_ode_functions(self,f,L):
        if not self.secondary_parent:
            ODE = sp.diff(sp.diff(L,self.dq1),self.t) - sp.diff(L,self.q1(self.t)) 
            f.append(ODE)
            ODE = sp.diff(sp.diff(L,self.dq2),self.t) - sp.diff(L,self.q2(self.t))
            f.append(ODE)

    def substitude_symbols(self, f):
        if not self.secondary_parent:
            for i in range(len(f)):
                f[i] = f[i].subs(sp.diff(self.q1(self.t),self.t,2),self.ddQ1).subs(sp.diff(self.q1(self.t),self.t),self.dQ1).subs(self.q1(self.t),self.Q1)
                f[i] = f[i].subs(sp.diff(self.q2(self.t),self.t,2),self.ddQ2).subs(sp.diff(self.q2(self.t),self.t),self.dQ2).subs(self.q2(self.t),self.Q2)
        

    def get_symbol(self):
        s=[]
        if not self.secondary_parent:
            s.append((self.Q1, self.dQ1, self.ddQ1))
            s.append((self.Q2, self.dQ2, self.ddQ2))
        return s

    def get_x0(self, x0):
        if not self.secondary_parent:
            x0.extend([self.x0, self.dx0])
            x0.extend([self.phi0, self.dphi0])

    def update(self, x):
        if not self.secondary_parent:
            self.x = x[2*self.x_index]
            self.dx = x[2*self.x_index+1]
            self.phi = x[2*self.phi_index]
            self.dphi = x[2*self.phi_index+1]

    def potential_expr(self,g):
        if not self.secondary_parent:
            return sp.Rational(1,2)*self.k*self.q1(self.t)**2
        pos1 = self.parent.get_position_expr()
        pos2 = self.secondary_parent.get_position_expr()
        x = sp.sqrt((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2)
        return sp.Rational(1,2)*self.k*x**2


    def plot(self, t=0):
        x,y = np.array([self.l2g(0,t),self.l2g(1,t)]).T
        self.plt_data.set_data(x,y)
        return self.plt_data

    def init_plot(self,ax):
        self.ax = ax
        self.plt_data, = ax.plot([],[],'-y')

    def evaluate(self,L):
        if self.secondary_parent:
            return L
        L = L.subs(self.dQ1,self.dx).subs(self.dq1,self.dx).subs(self.Q1,self.x).subs(self.q1(self.t),self.x)
        return L.subs(self.dQ2,self.dphi).subs(self.dq2,self.dphi).subs(self.Q2,self.phi).subs(self.q2(self.t),self.phi)