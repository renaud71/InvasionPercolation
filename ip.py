import numpy as np
from scipy.ndimage import measurements
import os

class IP:

	def __init__(self,Ncol,Nrow,theta):

		self.Ncol = Ncol
		self.Nrow = Nrow
		self.theta = theta
		self.dimensional_constant = 0.140/self.Nrow
		self.g = 9.81
		self.rho_nw = 1.225
		self.rho_w = 1205
		self.L = 1
		self.path='/Users/Toussaint/Desktop/BureauFracflow/Echanges/Monem/Invasion-Percolation-Python-master/IP/'
		self.pathsimulation = self.path+"Ncol"+str(self.Ncol)+"_Nrow"+str(self.Nrow)+"_theta"+str(self.theta)+"/"
		self.pmin = 200
		self.pmax = 608
		if not os.path.exists(self.pathsimulation):
			os.makedirs(self.pathsimulation)

		self.f = open(self.pathsimulation+ "widthofthefront.txt",'w') #

	def structure(self):
		"""
		Constructing the porous matrix
		"""
		L,g,Ncol,Nrow = self.L,self.g,self.Ncol,self.Nrow
		theta = (self.theta*np.pi)/180
		#Making the random numbers predictable
		np.random.seed(1)

		#initializing pressure and invasion lattice
		pressurelattice = np.zeros([Nrow,Ncol],float)
		invasionlattice = np.ones([Nrow,Ncol],int)
		delta_rho = self.rho_w-self.rho_nw


		#filling the lattice.
		for i in xrange(Nrow):
			for j in xrange(Ncol):
				pressurelattice[i,j]= np.random.uniform(self.pmin,self.pmax) - self.dimensional_constant*(Nrow-i)*g*delta_rho*np.sin(theta)
		

		return pressurelattice,invasionlattice

	def invasion(self):

		Ncol,Nrow = self.Ncol,self.Nrow

		#setting up the lattice
		pressurelattice, invasionlattice = self.structure()
		pressure = []
		number_of_clusters = []

		#boundary condition
		invasionlattice[0,:] = 0

		#Temporal plotting
		self.temporalplot = np.zeros([Nrow,Ncol],int)
		count = 1; tmp = self.pmax+100; xtmp = 0; ytmp = 0;

		#Time loop
		while True:

			#looks for trapped clusters hack.
			copy = np.copy(invasionlattice[-1,:])
			invasionlattice[-1,:] = 1
			lw, num = measurements.label(invasionlattice)
			number_of_clusters.append(num)

			lw[-1,:] = copy*lw[-1,:]
			invasionlattice[-1,:] = copy
			self.lw = lw

			#looping through the lattice.
			for i in xrange(Nrow):
				sumrow = np.sum(invasionlattice[i,:])
				if sumrow < Ncol and i != (self.Nrow-1):
					for j in xrange(Ncol):
						boolean = self.interfacesite(y=i+1,x=j)
						if boolean and pressurelattice[i+1,j]<tmp: #her ligger problemet
							tmp  = pressurelattice[i+1,j]; ytmp = i+1; xtmp = j


			invasionlattice[ytmp,xtmp] = 0
			pressure.append(tmp)
			number_of_clusters.append(num)
			self.frontwidth(invasionlattice)

			tmp = self.pmax+100;
			self.temporalplot = self.temporalplot + invasionlattice	
			count += 1;

			#Cuts off the simulation once the filter has been saturated
			if np.max(invasionlattice[-1,:]) == 0:
				print "Number of timesteps= %s"%(count)
				break

		lw, num = measurements.label(invasionlattice)
		#write to file
		np.savetxt(self.pathsimulation + "pressure"+str(Ncol)+"_"+str(Nrow)+".txt",[pressure])
		np.savetxt(self.pathsimulation + "bw"+str(Ncol)+"_"+str(Nrow)+".txt",invasionlattice)
		np.savetxt(self.pathsimulation + "temporal"+str(Ncol)+"_"+str(Nrow)+".txt",self.temporalplot)
		np.savetxt(self.pathsimulation + "clustersize"+str(Ncol)+"_"+str(Nrow)+".txt",lw)

		#plotting figures
		self.plotting(pressurelattice,invasionlattice,pressure,number_of_clusters)

	def frontwidth(self,invasionlattice):

		#This function is suppose to calculate the width of the front.
		for i in xrange(self.Ncol):
			for j in xrange((self.Nrow-1),-1,-1):
				if invasionlattice[j,i] == 0:
					self.f.write(str(j)+" ")
					break

		self.f.write("\n")

	def interfacesite(self,y,x):

		lw = self.lw
		value = lw[y,x]
		lv = 0
		#excluding trapped sites 
		if value != max(lw[-1,:]):
			return False
			
			#inner points
		elif  0 < x <(self.Ncol-1) and 0< y <(self.Nrow-1):
			if (lw[y-1,x] == lv or lw[y+1,x] == lv or lw[y,x-1]==lv or lw[y,x+1] ==lv):
				return True
			else:
				return False

			#Right boundary	
		elif x == (self.Ncol-1):
			if y != (self.Nrow-1) and (lw[y+1,x]==lv or lw[y-1,x]==lv or lw[y,x-1]==lv):
				return True
			elif y == (self.Nrow-1) and (lw[y-1,x]==lv or lw[y,x-1]==lv):
				return True 
			else:
				return False

			#Left boundary
		elif x == 0:
			if y != (self.Nrow-1) and (lw[y+1,x] ==lv or lw[y-1,x]==lv or lw[y,x+1]==lv):
				return True
			if y== (self.Nrow-1) and  (lw[y-1,x]==lv or lw[y,x+1]==lv):
				return True
			else:
				return False

			#Bottom boundary
		elif y == (self.Nrow-1):
			if (lw[y-1,x]==lv or lw[y,x-1]==lv or lw[y,x+1]==lv):
				return True
			#disregarded points
			else:
				return False
		else:
			return False


	def plotting(self,pressurelattice,invasionlattice,pressure,number_of_clusters):

		from matplotlib.pylab import plt

		#Plotting the invasionlattice
		plt.figure(2)
		plt.imshow(invasionlattice, cmap='Greys',interpolation='nearest')
		plt.title("Invasion lattice")
		plt.colorbar()
		plt.savefig(self.pathsimulation + "invasionlattice.png", bbox_inches="tight")
		plt.close()

		#plotting the pressure
		plt.figure(5)
		plt.plot(pressure)
		plt.xlabel('Time')
		plt.ylabel('Pressure')
		plt.title('P(t)')
		plt.savefig(self.pathsimulation +"pressure.png", bbox_inches="tight")
		plt.close()

		#plotting the clusterlattice
		plt.figure(6)
		plt.imshow(self.lw, interpolation='nearest')
		plt.title("Clusterlattice")
		plt.colorbar()
		plt.savefig(self.pathsimulation +"clusterlattice.png", bbox_inches="tight")
		plt.close()

		#Temporal diagram
		plt.figure(7)
		plt.imshow(self.temporalplot,interpolation='nearest')
		plt.title('Temporal diagram')
		plt.colorbar()
		plt.savefig(self.pathsimulation +"temporal.png", bbox_inches="tight")
		plt.close()

		#Plotting pressure distribution in the cell.
		plt.figure(3)
		plt.hist(pressurelattice.ravel(), bins=30, fc='k', ec='k')
		plt.savefig(self.pathsimulation +"pressuredistribution.png", bbox_inches="tight")
		plt.close()

		#Plotting the number of clusters as a function of interation.
		plt.figure(1)
		plt.plot(number_of_clusters)
		plt.xlabel('Iteration number/time')
		plt.ylabel('Number of clusters')
		plt.title('Number_of_cluster(t)')
		plt.savefig(self.pathsimulation +"clusterN.png", bbox_inches="tight")
		plt.close()

def execute():

	#Checking the execution time.
	import time
	start_time = time.time()
	c=IP(Nrow=7, Ncol=5, theta=15)
	c.invasion()
	print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == '__main__':
	execute()
