import numpy as np

Nx = 120; Ny=240; theta=0
a =  "Nx"+str(Nx)+"_"+"Ny"+str(Ny)+"_theta"+str(theta)
path = '/home/guillaume/Documents/python/IP/'+a+'/'
b = str(Nx)+"_"+str(Ny)
#The pressure 
pressure = np.loadtxt(path+'pressure'+b+'.txt')
Nsites = float(Nx*Ny)
steps = np.arange(1,(len(pressure)+1),1)
saturation = steps/Nsites
#Binary image
bw = np.loadtxt(path +'bw'+b+'.txt')
#clusters
clusters =  np.loadtxt(path +'clustersize'+b+'.txt')
#temporal
temporal = np.loadtxt(path + 'temporal'+b+'.txt')


#plotting.
import matplotlib.pylab as plt
plt.figure(1)
plt.plot(saturation,pressure,'.-')
plt.xlabel('saturation')
plt.ylabel('Pressure')

plt.figure(2)
plt.imshow(bw, cmap='Greys',interpolation='nearest')
plt.title("Invasion lattice")
plt.colorbar()

plt.figure(3)
plt.imshow(clusters,interpolation='nearest')
plt.title("Clusters")
plt.colorbar()

plt.figure(4)
plt.imshow(temporal,interpolation='nearest')
plt.title("Temporal evolution")
plt.colorbar()
plt.show()
