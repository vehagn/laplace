from __future__ import division
from numpy import *
from matplotlib import *
from pylab import *
from scipy import *
from scipy.integrate import quad

def stepfunc(x):
	if x > 0:
		return 1.
	else:
		return 0.
		
def deltapot(x):
	if x > 0.75 and x < 0.788:
		return 1.
	else:
		return 0.

def compute(pot=1,terms=20,tiles=100,potperiod=1,clines=30,fig=1,filename='fig',ext='.png'):
	
	#pot choses the true potental which from the ones given in V0
	#terms is the number of terms which is wanted in the approximation, it is not advisable to put this to more than 200.
	#tiles is the dimension of the grid which is used to calculate potential, increasing this increases the resolution, must be a multiple of 20 and minimum 40.
	#potperiod is the number of period for the sine-potential, only usable with pot = 1
	#clines are the number of contourlines which is to be drawn on the graphical representation of the potential
	#fig determines the window name in which the figure is drawn, must be a number.
	#filename is the name of the file in which the graph is saved.
	#ext is the extension og the file in ehich the graph is saved. Common extensions include .png and .pdf.
	
	if (tiles%20 != 0 ) | (tiles < 40):
		return '"tiles" must be a multiple of 20 and minimum 40'

	clf()
	C = range(terms)		
	x = y = linspace(0,1,tiles)
	X,Y = meshgrid(x,y)
	grid = eye(tiles)
	
	def V0(x):
		if pot == 1:
			return sin(2*potperiod*pi*x)
		if pot == 2:
			return 1-(x-1/2)**4
		if pot == 3:
			return stepfunc(x-1/2)*stepfunc(3/4-x)
		if pot == 4:
			return 2*sin(5.5*pi*x+1) + 3*stepfunc(x-1/8)*stepfunc(2/8-x) - 4*stepfunc(x-3/8)*stepfunc(4/8-x) - 2*stepfunc(x-3/8)*stepfunc(5/8-x) + 5*stepfunc(x-6/8)*stepfunc(7/8-x)
		if pot == 5:
			return 1/sqrt(2)*(sin(1*pi*x)+sin(4*pi*x))
		else:
			return x**0;
		
	for n in range(terms):
		(I,err) = quad(lambda x:V0(x)*sin((n+1)*pi*x),0,1,limit=100)
		C[n] = (2/sinh((n+1)*pi))*I #The Fourier-coefficients.
	
	def V(x,y):
		#Generates the model potential.
		Vm = 0;
		for n in range(terms):
				Vm += C[n]*sin((n+1)*pi*x)*sinh((n+1)*pi*y)
		return Vm
			
	grid = V(X,Y)
	
	figure(fig)
	setp(axes([0.20, 0.16, 0.64, 0.64]),xticks=[],yticks=[])
	contour(X,Y,grid,clines)
	xlabel('x/L')
	ylabel('y/L')

	xcomp = eye(19)	#Becomes the x-component of the E-field.
	ycomp = eye(19)	#Becomes the y-component of the E-field.
	
	for i in range(int(tiles/20),tiles,int(0.05*tiles)):
		for j in range(int(tiles/20),tiles,int(0.05*tiles)):
			a,b = int(20*i/tiles-1),int(20*j/tiles-1)
			xcomp[a,b] = (grid[i,j-1] - grid[i,j+1])
			ycomp[a,b] = (grid[i-1,j] - grid[i+1,j])
			
	chi = [x[i] for i in range(int(tiles/20),tiles,int(0.05*tiles))]
	psi = [y[j] for j in range(int(tiles/20),tiles,int(0.05*tiles))]
	Chi, Psi = meshgrid(chi,psi)
	quiver(Chi,Psi,xcomp,ycomp)
	
	V0_graf = range(len(x))
	SSE = 0;
	
	for i in range(len(x)):	#Stupid way to do it, but got errors otherwise for the stepfunc-dependent potentials.
		V0_graf[i] = V0(x[i])
		SSE += sqrt((V0_graf[i]-V(x[i],1))**2)
	
	setp(axes([0.20, 0.82, 0.64, 0.09]), xticks=[], yticks=[])
	title('V(x,L)')
	plot(x,V0_graf,'r',label='true')
	plot(x,V(x,y[-1]),'b',label='model')
	legend(ncol=2,loc=(.0,1.0),frameon=0,columnspacing=0.2,handletextpad=0.1)
	axhspan(0,0)

	setp(axes([0.12, 0.16, 0.04, 0.64]),xticks=[])
	xlim(-2e-15, 2e-15)
	xlabel('$\pm$ 2e-15')
	title('V(0,y)')
	plot(x-x,y,'r')
	plot(V(0,y),y,'b')
	
	setp(axes([0.88, 0.16, 0.04, 0.64]),xticks=[])
	xlim(-2e-15, 2e-15)
	xlabel('$\pm$ 2e-15')
	title('V(L,y)')
	plot(x-x,y,'r')
	plot(V(x[-1],y),y,'b')
	
	setp(axes([0.20, 0.08, 0.64, 0.04]),yticks=[])
	ylim(-3e-16,3e-16)
	xlabel('V(x,0)')
	plot(x,y-y,'r')
	plot(x,V(x,0),'b')
	
	savefig(filename+ext)
	
	print 'max abs V(L,y):    ', max(abs(V(x[-1],y)))
	print 'max abs V(x[-2],y): ', max(abs(V(x[-2],y)))	
	
	if SSE/tiles < 1e-10:
		print	'\nNo convergence issues'
	elif SSE/tiles <	0.05:
		print	'\nLight convergence issues, consider adding more terms'
	elif SSE/tiles <	0.1:
		print	'\nModerate convergence issues, consider adding more terms'
	else:
		print 	'\nSevere convergence issues, consider adding more terms'
