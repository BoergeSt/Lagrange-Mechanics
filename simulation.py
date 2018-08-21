
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import scipy.integrate
import numpy as np
import sympy as sp

from time import time

from sympy.utilities.iterables import flatten




class Simulation:
    def __init__(self,dt = 1./30,movie=False,subintegrations = 10):
        self.dt = dt
        self.Objects = []
        self.fig = plt.figure(figsize=(20,10))
        self.xlim = (-2,2)
        self.ylim = (-2,2)
        self.subintegrations = subintegrations
        
        self.ax = self.fig.add_subplot(111,aspect='equal', autoscale_on=False, xlim=self.xlim, ylim=self.ylim)
        self.movie = movie


    def addObjects(self,objects):
        self.Objects.extend(objects)

    def plot(self):
        self.ax.cla()
        self.ax.autoscale(False)
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)

        for object in self.Objects:
            object.plot(self.ax);

    def setup(self):
        i = 0;
        t = sp.symbols('t')
        for object in self.Objects:
            i = object.setup(i,t)
            

    
    def calculate_potential_expr(self):
        U = 0;
        for i,object in enumerate(self.Objects):
            Ui = object.potential_expr() 
            print("Object {} has potential {}".format(i,Ui))
            U += Ui

        return sp.simplify(U)

    def calculate_kinetic_expr(self):
        T = 0;
        for i,object in enumerate(self.Objects):
            Ti = object.kinetic_expr()
            print("Object {} has kinetic Energy {}".format(i,Ti))
            T+=Ti
        return sp.simplify(T)

    def calculate_lagrange_expr(self):
        L = sp.simplify(self.calculate_kinetic_expr()-self.calculate_potential_expr())
        print("Lagrange Function: L = {}".format(L))
        return L

    def calculate_ode_functions(self,L):
        f = []
        for object in self.Objects:
            #print("befor: {}".format(f))
            f = object.calculate_ode_functions(f,L)
            #print("after: {}".format(f))
        for object in self.Objects:
            f = object.substitude_symbols(f)

        s = self.get_symbols()
        for i in range(len(s)):
            s[i]=s[i][2]
        #print(s)
    
        rf = sp.solve(f,s)
        #print(rf[s[0]].free_symbols)
        for i in range(len(f)):
            f[i]=rf[s[i]]
        return f;

    def get_rhs(self,f,x):
        res = []
        for i in range(len(f)):
            res.extend([x[2*i+1],f[i]])
            for object in self.Objects:
                res[2*i+1] = object.substitude_values(res[2*i+1],x);
            #print("rhs[{}]={}\nrhs[{}]={}".format(2*i,res[2*i],2*i+1,res[2*i+1]))
        return res

    def get_x0(self):
        x0 = []
        for object in self.Objects:
            object.get_x0(x0);
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
        self.plot()
    
        self.setup()
        L = self.calculate_lagrange_expr()
        f = self.calculate_ode_functions(L)
        #print(f)
        #print(sim.get_symbols())
    
        s=self.get_symbols()
        for i in range(len(s)):
            s[i]=(s[i][0],s[i][1])

        #fl = []
        #for i,fi in enumerate(f):
        #	fl.append(sp.lambdify(s,s[i][0]))
        #	fl.append(sp.lambdify(s,fi))
        #print(fl)
        f2 = []
        for i,fi in enumerate(f):
            f2.extend([s[i][1],fi])
        print("ODE System: {}".format(f2))
        func = sp.lambdify(flatten(s),f2)
        rhs2 = lambda t,x: func(*x)
    


    
        #rhs = lambda t,x: sim.get_rhs(f,x)
        x0 = self.get_x0()
        #print(x0)
    
        #t0 = time()
        #rhs(0,x0)	
        #t1 = time()
        #print("rhs took {} ms".format(t1-t0))
        #t0 = time()
        #rhs2(0,x0)	
        #t1 = time()
        #print("rhs2 took {} ms".format(t1-t0))
    
        r = scipy.integrate.ode(rhs2).set_integrator('vode', method='bdf',with_jacobian=False) #bdf/adams
        r.set_initial_value(x0,0)


        def animate(i):
            for i in range(self.subintegrations):
                r.integrate(r.t+self.dt/self.subintegrations)
            self.update(r.y)
            #print(r.t)
            self.plot()


    
        t0 = time()
        animate(0)
        t1 = time()
        interval = 1000 * self.dt - (t1 - t0)


        anim = animation.FuncAnimation(self.fig, animate, frames=300,
                                      interval=interval)

        if self.movie:
            anim.save('lagrange.mp4', fps=30, extra_args=['-vcodec', 'libx264'])


        plt.show()
        
