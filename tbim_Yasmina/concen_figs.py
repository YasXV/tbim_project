import numpy as np
import matplotlib.pyplot as plt 

concen6=np.loadtxt("concentration6.dat",dtype=float)#Loading the concentration of type 6 
concen7=np.loadtxt("concentration7.dat",dtype=float)#Loading the concentration of type 7 
concen8=np.loadtxt("concentration8.dat",dtype=float)#Loading the concentration of type 8 
concen9=np.loadtxt("concentration9.dat",dtype=float)#Loading the concentration of type 9 
concen12=np.loadtxt("concentration12.dat",dtype=float)#Loading the concentration of type 12 

plt.scatter(concen6[:,0],concen6[:,1],label='concentration6',color='red')
plt.scatter(concen7[:,0],concen7[:,1],label='concentration7',color='green')
plt.scatter(concen8[:,0],concen8[:,1],label='concentration8',color='blue')
plt.scatter(concen9[:,0],concen9[:,1],label='concentration9',color='deeppink')
plt.scatter(concen12[:,0],concen12[:,1],label='concentration12',color='orange')
plt.legend()
plt.xlabel('dmu')
plt.ylabel('concentration')
plt.grid()


