
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import scipy.integrate
import numpy as np
import sympy as sp

from time import time

from sympy.utilities.iterables import flatten

import logging

logger = logging.getLogger('Lagrange_Mechanics')
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler('mechanics.log')
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
logger.addHandler(fh)
logger.addHandler(ch)


class Simulation:
    def __init__(self,dt = 1./30,movie=False,subintegrations = 10,xlim = (-2,2),ylim=(-2,2), g = scipy.constants.g, show_information = True):
        self.dt = dt
        self.Objects = []
        self.fig = plt.figure(figsize=(20,10))
        self.xlim = xlim
        self.ylim = ylim
        self.subintegrations = subintegrations
        self.t = sp.symbols('t')
        
        self.ax = self.fig.add_subplot(111,aspect='equal', autoscale_on=False, xlim=self.xlim, ylim=self.ylim)
        self.movie = movie
        logger.debug("Simulation initiated.")
        self.g = g
        self.show_information = show_information


    def addObjects(self,objects):
        self.Objects.extend(objects)
        logger.debug("Added objects")

    def init_plot(self):
        for object in self.Objects:
            object.init_plot(self.ax)

    def plot(self, t = 0):
        #self.ax.cla()
        #self.ax.autoscale(False)
        #self.ax.set_xlim(self.xlim)
        #self.ax.set_ylim(self.ylim)
        res = []

        for object in self.Objects:
            res.append(object.plot(t))
        return res

    def setup(self): #TODO: can fail if objects are in wrong order
        logger.debug("Starting setup")
        i = 0 
        for object in self.Objects:
            i = object.setup(i,self.t)
        logger.debug("Finished setup of {} objects with {} independent variables".format(len(self.Objects),i))
            

    
    def calculate_potential_expr(self):
        U = 0 
        for i,object in enumerate(self.Objects):
            Ui = object.potential_expr(self.g) 
            logger.debug("Object {} has potential {}".format(i,Ui))
            U += Ui

        return sp.simplify(U)

    def calculate_kinetic_expr(self):
        T = 0 
        for i,object in enumerate(self.Objects):
            Ti = object.kinetic_expr()
            logger.debug("Object {} has kinetic Energy {}".format(i,Ti))
            T+=Ti
        return sp.simplify(T)

    def calculate_lagrange_expr(self):
        T = self.calculate_kinetic_expr()
        U = self.calculate_potential_expr()
        L = sp.simplify(T-U)
        H = sp.simplify(T+U)
        Lp = [L,H]
        for object in self.Objects:
            object.substitude_symbols(Lp)
        logger.info("Lagrange Function:\n\tL = {}".format(Lp[0]))
        logger.info("Lagrange Function unchanged:\n\tL = {}".format(L))
        return L,Lp[0],Lp[1]

    def calculate_ode_functions(self,L):
        f = []
        for object in self.Objects:
            object.calculate_ode_functions(f,L)
        for object in self.Objects:
            object.substitude_symbols(f)

        s = self.get_symbols()
        for i in range(len(s)):
            s[i]=s[i][2]
    
        logger.debug("cuppled ode: {}".format(f))
        rf = sp.solve(f,s)
        logger.debug("solver returned: {}".format(rf))
        for i in range(len(f)):
            try:
                f[i] = rf[s[i]]
            except KeyError:
                f[i] = sp.Integer(0)
                logger.warning("Symbol {} not found in solution. Setting it to 0".format(s[i]))
        for i,fun in enumerate(f):
            f[i] = sp.simplify(fun)
        return f 

    def get_x0(self):
        x0 = []
        for object in self.Objects:
            object.get_x0(x0) 
        return x0

    def get_symbols(self):
        s = []
        for object in self.Objects:
            s.extend(object.get_symbol())
        return s

    def update(self,x):
        for object in self.Objects:
            object.update(x)

    def evaluate(self,L,t):
        for object in self.Objects:
            L = object.evaluate(L)
        L = L.doit()
        L = L.subs(self.t,t)
        return L
    

    def run(self):
        self.setup()
        self.init_plot()
        self.plot()
        L,Lp,H = self.calculate_lagrange_expr()
        f = self.calculate_ode_functions(L)
    
        s=self.get_symbols()
        for i in range(len(s)):
            s[i]=(s[i][0],s[i][1])

        f2 = []
        for i,fi in enumerate(f):
            f2.extend([s[i][1],fi])
        logger.info("ODE System: \n\t{}".format("\n\t".join([str(o) for o in f2])))
        func = sp.lambdify([self.t]+flatten(s),f2)
        rhs2 = lambda t,x: func(t,*x)
    
        x0 = self.get_x0()
        logger.debug("x0 = {}".format(x0))

    
        r = scipy.integrate.ode(rhs2).set_integrator('vode', method='adams',with_jacobian=False) #bdf/adams
        r.set_initial_value(x0,0)

        logger.debug("Initialized Integrator")

        time_template = 'time = %.1fs'
        energy_template = 'energy = %.1f'
        time_text = self.ax.text(0.05, 0.9, '', transform=self.ax.transAxes)
        energy_text = self.ax.text(0.05,0.87,'', transform=self.ax.transAxes)


        def animate(i):
            for i in range(self.subintegrations):
                r.integrate(r.t+self.dt/self.subintegrations)
            if not r.successful():
                r.t += self.dt
            self.update(r.y)
            res = self.plot(r.t)
            if self.show_information:
                time_text.set_text(time_template % r.t)
                energy_text.set_text(energy_template % (self.evaluate(H,r.t)))
                res.extend([time_text,energy_text])
            return tuple(res)


    
        t0 = time()
        animate(0)
        t1 = time()
        interval = 1000 * self.dt - (t1 - t0)


        anim = animation.FuncAnimation(self.fig, animate, frames=300,
                                      interval=interval, blit = True)

        if self.movie:
            anim.save('lagrange.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

        logger.debug("Starting animation")
        plt.show()
        logger.debug("Finished animation")
        
