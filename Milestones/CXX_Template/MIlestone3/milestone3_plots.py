import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("../perturbations_k0.01.txt")
#x, theta0, theta1, Phi, Psi, vb, v_cdm, delta_b, delta_cdm
"""
n = int(len(data[:,0]))
x = data[:,0]
z = 1/np.exp(x) - 1

theta0 = data[:,1]
theta1 = data[:,2]

Phi = data[:,3]
Psi = data[:,4]

vb = data[:,5]
v_cdm = data[:,6]
delta_b = data[:,7]
delta_cdm = data[:,8]

plt.plot(x,theta0)
plt.title("theta0")
plt.show()

plt.plot(x,theta1)
plt.title("theta1")
plt.show()

plt.plot(x,Phi)
plt.title("phi")
plt.show()

plt.plot(x,Psi)
plt.title("psi")
plt.show()

plt.plot(x,vb)
plt.plot(x,v_cdm)
plt.title("vb,vcdm")
plt.show()

plt.plot(x,delta_b)
plt.plot(x,delta_cdm)
plt.title("deltas")
plt.show()
"""

