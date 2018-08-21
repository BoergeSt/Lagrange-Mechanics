
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
    def __init__(self,dt = 1./30,movie=False,subintegrations = 10,xlim = (-2,2),ylim=(-2,2)):
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


    def addObjects(self,objects):
        self.Objects.extend(objects)
        logger.debug("Added objects")

    def plot(self, t = 0):
        self.ax.cla()
        self.ax.autoscale(False)
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)

        for object in self.Objects:
            object.plot(self.ax,t) 

    def setup(self):
        logger.debug("Starting setup")
        i = 0 
        for object in self.Objects:
            i = object.setup(i,self.t)
        logger.debug("Finished setup of {} objects with {} independent variables".format(len(self.Objects),i))
            

    
    def calculate_potential_expr(self):
        U = 0 
        for i,object in enumerate(self.Objects):
            Ui = object.potential_expr() 
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
        L = sp.simplify(self.calculate_kinetic_expr()-self.calculate_potential_expr())
        logger.info("Lagrange Function:\n\tL = {}".format(L))
        return L

    def calculate_ode_functions(self,L):
        f = []
        for object in self.Objects:
            object.calculate_ode_functions(f,L)
        for object in self.Objects:
            object.substitude_symbols(f)

        s = self.get_symbols()
        for i in range(len(s)):
            s[i]=s[i][2]
    
        rf = sp.solve(f,s)
        for i in range(len(f)):
            f[i]=rf[s[i]]
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

    def run(self):
        self.setup()
        self.plot()
        L = self.calculate_lagrange_expr()
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

    
        r = scipy.integrate.ode(rhs2).set_integrator('vode', method='bdf',with_jacobian=False) #bdf/adams
        r.set_initial_value(x0,0)

        logger.debug("Initialized Integrator")


        def animate(i):
            for i in range(self.subintegrations):
                r.integrate(r.t+self.dt/self.subintegrations)
            self.update(r.y)
            #print(r.t)
            self.plot(r.t)


    
        t0 = time()
        animate(0)
        t1 = time()
        interval = 1000 * self.dt - (t1 - t0)


        anim = animation.FuncAnimation(self.fig, animate, frames=300,
                                      interval=interval)

        if self.movie:
            anim.save('lagrange.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

        logger.debug("Starting animation")
        plt.show()
        logger.debug("Finished animation")
        
