import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 14})
import os

#print(os.listdir("../"))



#DATA FROM PERTURBATION CLASS
data_03 = np.genfromtxt("../perturbations_k0.3.txt")
data_01 = np.genfromtxt("../perturbations_k0.1.txt")
data_001 = np.genfromtxt("../perturbations_k0.01.txt")
data_0001 = np.genfromtxt("../perturbations_k0.001.txt")
data_00005 = np.genfromtxt("../perturbations_k0.0005.txt")
#x, theta0, theta1, theta2,  Phi, Psi, vb, v_cdm, delta_b, delta_cdm, source_function
data_list = [data_03, data_01, data_001, data_0001,data_00005]
k_list = [0.3,0.1,0.01,0.001,0.0005]
colors = ["royalblue","darkorange","seagreen", "darkorchid", "maroon"]

n = int(len(data_01[:,0]))
x = data_01[:,0]

for i in range(len(data_list)):
    source_function_i = data_list[i][:,10]
    plt.plot(x, source_function_i, label = f"k={k_list[i]}")
plt.legend()
plt.show()