import numpy as np
import matplotlib.pyplot as plt


data_01 = np.genfromtxt("../perturbations_k0.1.txt")
data_001 = np.genfromtxt("../perturbations_k0.01.txt")
data_0001 = np.genfromtxt("../perturbations_k0.001.txt")
data_list = [data_01, data_001, data_0001]
#x, theta0, theta1, Phi, Psi, vb, v_cdm, delta_b, delta_cdm

n = int(len(data_01[:,0]))
x = data_01[:,0]
z = 1/np.exp(x) - 1

def plot_theta0(x, data_list):
    for data in data_list:
        theta0 = data[:,1]
        plt.plot(x, theta0)
        plt.title("Theta0")
    plt.legend(["k=0.1", "k=0.01", "k=0.001"])
    plt.show()

plot_theta0(x, data_list)


def plot_theta1(x, data_list):
    for data in data_list:
        theta1 = data[:,2]
        plt.plot(x, theta1)
        plt.title("Theta1")
    plt.legend(["k=0.1", "k=0.01", "k=0.001"])
    plt.show()

plot_theta1(x, data_list)

def plot_Phi(x, data_list):
    for data in data_list:
        phi = data[:,3]
        plt.plot(x, phi)
        plt.title("Phi")
    plt.legend(["k=0.1", "k=0.01", "k=0.001"])
    plt.show()

plot_Phi(x, data_list)

def plot_Psi(x, data_list):
    for data in data_list:
        psi = data[:,4]
        plt.plot(x, psi)
        plt.title("Psi")
    plt.legend(["k=0.1", "k=0.01", "k=0.001"])
    plt.show()

plot_Psi(x, data_list)

def plot_vb_v_cdm(x, data_list):
    for data in data_list:
        vb = data[:,5]
        vb_cdm = data[:,6]
        plt.plot(x, vb,"--")
        plt.plot(x,vb_cdm)
        plt.title("v")
    plt.legend(["k=0.1", "k=0.01", "k=0.001"])
    plt.show()

plot_vb_v_cdm(x, data_list)

def plot_deltas(x, data_list):
    for data in data_list:
        deltab = data[:,7]
        delta_cdm = data[:,8]
        plt.plot(x, deltab,"--")
        plt.plot(x,delta_cdm)
        plt.title("delta")
    plt.legend(["k=0.1", "k=0.01", "k=0.001"])
    plt.show()
        

plot_deltas(x, data_list)
"""
theta0 = data[:,1]
theta1 = data[:,2]
"""

   

"""
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


