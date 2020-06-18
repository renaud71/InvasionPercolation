
import numpy as np
from scipy.ndimage import measurements
import os
import time

class IPfilm:

	def __init__(self,Nrows,Ncol,theta):

		self.Nrows = int(Nrows) #row
		self.Ncol = int(Ncol) #col
		self.theta = theta # convertion from angles to radians is needed.
		self.g = 1
		self.rho = 1
		self.L = 1
		self.path='/home/guillaume/Documents/python/IPfilm/'
		self.pathsimulation = self.path+"Nrows"+str(self.Nrows)+"_Ncol"+str(self.Ncol)+"_theta"+str(self.theta)+"/"

		if not os.path.exists(self.pathsimulation):
			os.makedirs(self.pathsimulation)
	def structure(self):

		L,g,Nrows,Ncol,rho,theta = self.L,self.g,self.Nrows,self.Ncol,self.rho, self.theta
		#element spcification.
		"""
		0 = liquid pore
		1 = bond
		2 = bead
		3 = invaded pore
		4 = film
		"""
		BL = np.zeros([Nrows+4,(Ncol+4)],int)
		PL = np.zeros([Nrows+4,Ncol+4],float)
		FILM = np.zeros([Nrows+4,(Ncol+4)],int)
		#fyller de ytre punktene med verdien 9.
		BL[:,-2:] = BL[-2:,:] = BL[:2,:] = BL[:,:2] = 9;
		#Bondary condition.
		for i in range(2,Nrows+2,1):
			if i%2 == 1:
				BL[i,2:-2:2] = 2
				BL[i,3:-2:2] = 1
			else:
				BL[i,2:-2:2] = 1
				BL[i,3:-2:2] = 0

		if BL[-3,2] == 2:
			BL[-4,3:-2:2] = 3
		else:
			BL[-3,3:-2:2] = 3

		#Making the random numbers predictable
		np.random.seed(1)
		#filling the presssure lattice.
		for i in xrange(Nrows+4):
			for j in xrange(Ncol+4):
				PL[i,j]= float(np.random.randint(100,1000))/1000 - (Nrows-i)*g*rho*np.sin(theta)

		self.BL = BL

		return BL,PL,FILM
	def invasion(self):
		#starting the clock
		start_time = time.time()

		Nrows, Ncol = self.Nrows, self.Ncol

		#setting up the lattices and pressure list.
		BL,PL,FILM = self.structure()
		self.FILM = FILM
		pressure = []; count = 1

		#Setting up temporal lattice.
		self.temporalplot = np.zeros([int(Nrows/2),int(Ncol/2)],int)

		while True:
			Pmax=1; itmp=0; jtmp=0; ibond=0; jbond=0

			SL = self.bond2site(BL)
			hack = np.copy(SL)-3
			lw, num = measurements.label(hack)

			clustercheck = self.site2bond(lw)

			#Looping trough the all the bonds in the bond matrix, checking the vertical situation of the inner points.
			for i in range(2,(Nrows+2),1):
				for j in range(2,(Ncol+2),1):

					#checking downwards
					if (BL[i+1,j] == 1 and BL[i+2,j] ==3 and BL[i,j] == 0 and PL[i+1,j]<Pmax) and clustercheck[i,j]==clustercheck[2,3]:
						Pmax = PL[i+1,j]
						itmp = i; jtmp=j; ibond = i+1; jbond=j

					#checking towards the right
					if (BL[i,j+1] == 1 and BL[i,j+2] == 3 and BL[i,j] == 0 and PL[i,j+1]<Pmax) and clustercheck[i,j]==clustercheck[2,3]:
						Pmax = PL[i,j+1]
						itmp = i; jtmp=j; ibond = i; jbond=j+1

					#checking towards the left.
					if	(BL[i,j-1] == 1 and BL[i,j-2] ==3 and BL[i,j] == 0 and PL[i,j-1]<Pmax) and clustercheck[i,j]==clustercheck[2,3]:
						Pmax = PL[i,j-1]
						itmp = i; jtmp=j; ibond = i; jbond=j-1

			BL[itmp,jtmp] = 3
			self.BL = BL
			pressure.append(Pmax)
			self.temporalplot = self.temporalplot + SL
			self.evaluatefilm(itmp,jtmp,ibond,jbond)
			count += 1

			#cutoff
			if sum(SL[0,:]) > 0:
				break

		#self.filmthresh()

		#write to file
		np.savetxt(self.pathsimulation + "pressure"+str(Nrows)+"_"+str(Ncol)+".txt",[pressure])
		np.savetxt(self.pathsimulation + "SL"+str(Nrows)+"_"+str(Ncol)+".txt",SL)
		np.savetxt(self.pathsimulation + "BL"+str(Nrows)+"_"+str(Ncol)+".txt",BL)
		np.savetxt(self.pathsimulation + "film"+str(Nrows)+"_"+str(Ncol)+".txt",FILM)
		np.savetxt(self.pathsimulation + "temporal"+str(Nrows)+"_"+str(Ncol)+".txt",self.temporalplot)
		np.savetxt(self.pathsimulation + "film"+str(Nrows)+"_"+str(Ncol)+".txt",lw)


		return BL, SL,FILM
	def bond2site(self,BL):

		Nrows, Ncol = self.Nrows, self.Ncol
		SL = np.zeros([int(Nrows/2),int(Ncol/2)],int)
		count = 0
		for i in range(2,Nrows+1,2):			
			if i%2 != 1:
				SL[count,:] = BL[i,3:-2:2]
				count += 1
			pass

		return SL
	def site2bond(self,SL):

		Nrows,Ncol = self.Nrows,self.Ncol
		BL = self.BL
		copy = np.copy(self.BL)

		count = 0
		for i in range(1,Nrows+1,1):
			if i%2 != 1:
				copy[i,3:-2:2] = SL[count,:]
				count+=1
		return copy
	def evaluatefilm(self,itmp,jtmp,ibond,jbond):

		BL = self.BL
		if (jtmp-jbond) < 0:
			#invasion to the left
			if (BL[itmp,jtmp] == BL[itmp+2,jtmp]):
				#Left
				self.BL[itmp+1,jtmp]  = 4
				self.FILM[itmp+1,jtmp] = 2		#Film
				self.FILM[itmp+1,jtmp-1] = 4	#Bead
				self.FILM[itmp+1,jtmp+1] = 4	#Bead					

			elif (BL[itmp,jtmp] == BL[itmp-2,jtmp]):
				#Right
				self.BL[itmp-1,jtmp]  = 4
				self.FILM[itmp-1,jtmp] = 2		#Film
				self.FILM[itmp-1,jtmp+1] = 4	#Bead
				self.FILM[itmp-1,jtmp-1] = 4	#Bead

			elif (BL[itmp,jtmp] == BL[itmp,jtmp-2]):
				#forward ?
				self.BL[itmp,jtmp-1]  = 4
				self.FILM[itmp,jtmp-1] = 2		#Film
				self.FILM[itmp+1,jtmp-1] = 4	#Bead
				self.FILM[itmp-1,jtmp-1] = 4	#Bead

		elif (itmp-ibond) < 0:
			#invasion from the bottom.
			if (BL[itmp,jtmp] == BL[itmp,jtmp+2]):
				#right
				self.BL[itmp,jtmp+1] = 4
				self.FILM[itmp,jtmp+1] = 2		#Film
				self.FILM[itmp+1,jtmp+1] = 4	#Bead
				self.FILM[itmp-1,jtmp+1] = 4	#Bead

			elif (BL[itmp,jtmp] == BL[itmp,jtmp-2]):
				#left
				self.BL[itmp,jtmp-1] = 4
				self.FILM[itmp,jtmp-1] = 2		#Film
				self.FILM[itmp+1,jtmp-1] = 4	#Bead
				self.FILM[itmp-1,jtmp-1] = 4	#Bead				

			elif (BL[itmp,jtmp] == BL[itmp-2,jtmp]):
				#up
				self.BL[itmp-1,jtmp] = 4
				self.FILM[itmp-1,jtmp] = 2		#Film
				self.FILM[itmp-1,jtmp+1] = 4	#Bead
				self.FILM[itmp-1,jtmp-1] = 4	#Bead		

		elif (jtmp-jbond) > 0:
			#invasion to the right
			if (BL[itmp,jtmp] == BL[itmp+2,jtmp]):
				#up
				self.BL[itmp+1,jtmp] = 4
				self.FILM[itmp+1,jtmp] = 2		#Film
				self.FILM[itmp+1,jtmp+1] = 4	#Bead
				self.FILM[itmp+1,jtmp-1] = 4	#Bead

			elif (BL[itmp,jtmp] == BL[itmp,jtmp+2]):
				#right
				self.BL[itmp,jtmp+1]  = 4
				self.FILM[itmp,jtmp+1] = 2		#Film
				self.FILM[itmp-1,jtmp+1] = 4	#Bead
				self.FILM[itmp+1,jtmp+1] = 4	#Bead

			elif (BL[itmp,jtmp] == BL[itmp-2,jtmp]):
				#down
				self.BL[itmp-1,jtmp]  = 4
				self.FILM[itmp-1,jtmp] = 2		#Film
				self.FILM[itmp-1,jtmp-1] = 4	#Bead
				self.FILM[itmp,jtmp+1] = 4		#Bead
	def plotting(self,BL, SL,FILM):
		import matplotlib.pylab as plt
		plt.figure(1)
		plt.imshow(self.BL[2:-2,2:-2], interpolation='nearest')
		plt.title("BL")
		plt.colorbar()
		plt.savefig(self.pathsimulation +"bondlattice.eps", bbox_inches="tight")

		plt.figure(2)
		plt.imshow(SL,cmap='Greys', interpolation='nearest')
		plt.title("SL")
		plt.colorbar()
		plt.savefig(self.pathsimulation +"pixellattice.eps", bbox_inches="tight")

		plt.figure(3)
		plt.imshow(self.FILM[2:-2,2:-2],cmap='Greys', interpolation='nearest')
		plt.title("film")
		plt.colorbar()
		plt.savefig(self.pathsimulation +"filmlattice.eps", bbox_inches="tight")


		plt.figure(4)
		plt.imshow(self.temporalplot,interpolation='nearest')
		plt.title("Temporal plot")
		plt.colorbar()
		plt.savefig(self.pathsimulation +"temporalplot.eps", bbox_inches="tight")
	def filmthresh(self):
		Nrows,Ncol = self.Nrows,self.Ncol
		BL = self.BL
		FILM = self.FILM
		#Allocating memory for filmflow lattice.
		for i in range(2,(Nrows+2),1):
			for j in range(2,(Ncol+2),1):
				if BL[i,j] == 4 or BL[i,j] == 2:
					FILM[i,j] = 1;
				else:
					FILM[i,j] = 0;

def execute():
	#Checking the execution time.
	import time
	start_time = time.time()
	c=IPfilm(400,200,0) # Nrow and Ncol must be even numbers.
	BL, SL,film = c.invasion()
	c.plotting(BL,SL,film)
	print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
	execute()
