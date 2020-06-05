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

terms_01 = np.genfromtxt("../sourcetermsk01.txt")
terms_001 = np.genfromtxt("../sourcetermsk001.txt")
terms_0001 = np.genfromtxt("../sourcetermsk0001.txt")

terms_list = [terms_01, terms_001, terms_0001]
k_list_terms = [0.1, 0.01, 0.001]
n_terms = int(len(x))
x_terms = terms_01[:,0]
print(len(x))
print(len(terms_list[0][:,1]))

for i in range(1):
    A_i = terms_list[i][:,1]
    B_i = terms_list[i][:,2]
    C_i = terms_list[i][:,3]
    D_i = terms_list[i][:,4]
    plt.plot(x_terms, A_i, label= f"A, k={k_list[i]}")
    plt.plot(x_terms, B_i, label= f"B, k={k_list[i]}")
    plt.plot(x_terms, C_i, label= f"C, k={k_list[i]}")
    plt.plot(x_terms, D_i, label= f"D, k={k_list[i]}")
plt.legend()
plt.show()

