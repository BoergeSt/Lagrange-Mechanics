import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import sympy as sp
import scipy.constants
import scipy.integrate
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
		


class FixPoint:
	def __init__(self,position=np.array([0,0])):
		self.position = position

	def get_position(self):
		return self.position

	def get_position_expr(self):
		return self.position

	def plot(self,ax):
		#print("Fixpoint at {}".format(self.get_position()))
		x,y = self.get_position().T
		ax.plot(x,y,'ok')

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

	def substitude_values(self,vi,x):
		return vi;

	def get_symbol(self):
		return []

	def get_x0(self,x0):
		pass

	def update(self,x):
		pass

class Point:
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

	def calculate_ode_functions(self,f,L):
		#print("in Point: {}".format(f))
		return f

	def substitude_symbols(self,f):
		return f

	def substitude_values(self,vi,x):
		return vi;

	def get_symbol(self):
		return []
	
	def get_x0(self,x0):
		pass

	def update(self,x):
		pass

class Connector:
	def __init__(self,parent,length=1,offset = 0,phi0 = 0, dphi0 = 0):
		self.parent = parent
		self.length = length
		self.offset = offset
		self.phi0 = phi0
		self.dphi0 = dphi0
		self.phi = phi0
		self.dphi= dphi0


	def l2g(self,local): #local to global
		return self.parent.get_position()+self.length*(self.offset+local)*np.array([np.sin(self.phi),-np.cos(self.phi)])

	def l2g_expr(self,local): #local to global expression
		modifier = self.length*(self.offset+local)
		par_pos = self.parent.get_position_expr()
		#print("l2g_expr: {}".format([par_pos[0]+modifier*sp.sin(self.q(self.t)),par_pos[1]-modifier*sp.cos(self.q(self.t))]))
		return [par_pos[0]+modifier*sp.sin(self.q(self.t)),par_pos[1]-modifier*sp.cos(self.q(self.t))]
	
	def plot(self,ax):
		x,y = np.array([self.l2g(0),self.l2g(1)]).T
		ax.plot(x,y,'-k')
		#print("Connector from {} to {}".format(self.l2g(0),self.l2g(1)))

	def setup(self,i,t):
		self.index = i
		self.q=sp.Function('q{}'.format(i))
		self.dq = sp.diff(self.q(t),t)
		self.t = t
		self.Q,self.dQ,self.ddQ = sp.symbols('Q{0} dQ{0} ddQ{0}'.format(self.index))
		return i+1

	def potential_expr(self):
		return 0

	def kinetic_expr(self):
		return 0

	def calculate_ode_functions(self,f,L):
		ODE = sp.diff(sp.diff(L,self.dq),self.t) - sp.diff(L,self.q(self.t)) 
		return f+[ODE]

#		ODE = ODE.subs(sp.diff(self.q(self.t),self.t,2),self.ddQ).subs(sp.diff(self.q(self.t),self.t),self.dQ).subs(self.q(self.t),self.Q)
#		
#		#print(ODE)
#		solution = sp.solve(ODE,self.ddQ)[0]
#		#print(solution);
#		f.extend([solution])
#		return f;
		
	def substitude_values(self,vi,x):
		#print("befor subs from object {} with {}: {}".format(self.index,x,vi))
		#print("substitution {}. In the equations are the symbols: {}".format([self.Q,self.dQ],vi.free_symbols))
		res =  vi.subs(self.Q,x[2*self.index]).subs(self.dQ,x[2*self.index+1]);
		#print("after subs from object {}: {}".format(self.index,res))
		return res

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
		self.dphi = x[2*self.index+1]



def example():
	sim = Simulation(movie=True)
	
	Base = FixPoint()
	C1 = Connector(Base,phi0=np.pi/4)
	P1 = Point(C1)
	sim.addObjects([Base,C1,P1])

	C2 = Connector(Base,phi0=np.pi/4)
	sim.addObjects([C2])
	Pn = [Point(C2,local = l,mass=1/50) for l in np.linspace(0,1,50)]
	sim.addObjects(Pn)
	

	C3 = Connector(Base,phi0=np.pi/4)
	sim.addObjects([C3])
	Pn2 = [Point(C3,local = l,mass=1/100) for l in np.linspace(0,1,100)]
	sim.addObjects(Pn2)

	C4 = Connector(Base, phi0 = np.pi/4)
	PL0 = Point(C4,local = (1+np.sqrt(3/5))/2,mass=5/18)
	PL1 = Point(C4,local = 1/2,mass=4/9)
	PL2 = Point(C4,local = (1-np.sqrt(3/5))/2,mass=5/18)
	sim.addObjects([C4,PL0,PL1,PL2])

	#C2 = Connector(P1,0.5)
	#P2 = Point(C2)
	#sim.addObjects([C2,P2])

	#C3 = Connector(P2,1/4)
	#P3 = Point(C3)
	#sim.addObjects([C3,P3])


	sim.run()


if __name__=="__main__":
	example()
	#import profile
	#profile.run("example()")


