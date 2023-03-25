import numpy as np
import matplotlib.pyplot as plt 

concen6=np.loadtxt("concentration6.dat",dtype=float)#Loading the concentration of type 6 
concen7=np.loadtxt("concentration7.dat",dtype=float)#Loading the concentration of type 7 
concen8=np.loadtxt("concentration8.dat",dtype=float)#Loading the concentration of type 8 
concen9=np.loadtxt("concentration9.dat",dtype=float)#Loading the concentration of type 9 
concen12=np.loadtxt("concentration12.dat",dtype=float)#Loading the concentration of type 12 

plt.plot(concen6[:,0],concen6[:,1],'o-',label='concentration6',color='red')
plt.plot(concen7[:,0],concen7[:,1],'o-',label='concentration7',color='green')
plt.plot(concen8[:,0],concen8[:,1],'-o',label='concentration8',color='blue')
plt.plot(concen9[:,0],concen9[:,1],'o-',label='concentration9',color='deeppink')
plt.plot(concen12[:,0],concen12[:,1],'o-',label='concentration12',color='orange')
plt.legend()
plt.xlabel(r'd$\mu$')
plt.ylabel('concentration in Ag')
plt.grid()
plt.show()


