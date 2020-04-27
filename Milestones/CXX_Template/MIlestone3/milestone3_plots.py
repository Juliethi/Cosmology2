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
#x, theta0, theta1, theta2,  Phi, Psi, vb, v_cdm, delta_b, delta_cdm
data_list = [data_03, data_01, data_001, data_0001,data_00005]
k_list = [0.3,0.1,0.01,0.001,0.0005]
colors = ["royalblue","darkorange","seagreen", "darkorchid", "maroon"]

n = int(len(data_01[:,0]))
x = data_01[:,0]
z = 1/np.exp(x) - 1


#DATA FROM BACKGROUND COSMOLOGY
cosmo_data = np.genfromtxt("../cosmology.txt")
#x, eta(x), Hp(x), dHp_dx, omegaB, omegaCDM, omegaLambda, OmegaR, OmegaNu, OmegaK"
Mpc = 3.08567758e22
#eta = cosmo_data[:,1]/Mpc/1000. #Converting eta from m to Mpc
eta = cosmo_data[:,1] #m
Hp = cosmo_data[:,2]*Mpc/1000. #Converting Hp to km/s/Mpc
H = Hp/np.exp(x)

omegaB = cosmo_data[:,4]
omegaCDM = cosmo_data[:,5]
omegaLambda = cosmo_data[:,6]
omegaR = cosmo_data[:,7]

omegaSum = omegaB + omegaCDM + omegaLambda + omegaR
omegaMatter = omegaCDM + omegaB

#DATA FROM RECOMBINATION CLASS
data_recombination = np.genfromtxt("../recombination.txt")
#x, Xe(X), Xe_Saha_only, ne(X), tau,tau', tau'', g,g',g''
g_tilde = data_recombination[:,7]
g_tilde_max = np.argmax(g_tilde) #Index where the visibility function peaks, which corresponds to when tau=1
Xe_half = np.argmin(abs(data_recombination[:,1]-0.5)) #Index when recombination is halfway

def find_eras():
    rad_matter_dif = abs(omegaMatter - omegaR)
    rad_matter_dif = rad_matter_dif[:int(n/2)]
    rad_matter_eq = np.argmin(rad_matter_dif)
    #print("x value rad matter eq:", x[rad_matter_eq])
    #print("z value of rad matter eq:", z[rad_matter_eq])
    
    matter_de_dif = abs(omegaMatter - omegaLambda)
    matter_de_dif = matter_de_dif[int(n/2):]
    matter_de_eq = np.argmin(matter_de_dif) + int(n/2)
    #print("x value matter de eq:", x[matter_de_eq])
    #print("z value of matter de eq:", z[matter_de_eq])

    return rad_matter_eq, matter_de_eq

rad_matter_eq, matter_de_eq = find_eras()

def find_horizon_entry():
    horizon_entry = []
    for k in k_list:
        k = k/Mpc
        entry = np.argmin(abs(k-1./eta))
        horizon_entry.append(entry)
    return horizon_entry
horizon_entry = find_horizon_entry()


def plot_theta0(x, data_list):
    y_min=-1
    y_max = 1.2
    for i in range(len(data_list)):
        theta0 = data_list[i][:,1]
        plt.plot(x, theta0, color = colors[i], label= f"k={k_list[i]}")
        plt.vlines(x[horizon_entry[i]], color = colors[i], ymin=y_min, ymax=y_max, linestyle="dashdot")
    plt.title(r"$\Theta_0$")
    plt.xlabel("x")
    plt.xlim(-12,-0)
    plt.ylim(y_min,y_max)
    plt.vlines(x[g_tilde_max], ymin=y_min, ymax=y_max, linestyle="--", label=r"$x_{rec}$")
    #plt.vlines(x[Xe_half], ymin=-0.75, ymax=1, linestyle="--", label=r"x_{rec}")
    plt.axvspan(x[0],x[rad_matter_eq],color="peachpuff")
    plt.axvspan(x[rad_matter_eq], x[matter_de_eq], color="powderblue")
    plt.axvspan(x[matter_de_eq], x[-1], color="pink")
    plt.legend(loc='lower right')
    #plt.tight_layout()
    plt.savefig("theta0_plot.pdf")
    plt.show()

plot_theta0(x, data_list)


def plot_theta1(x, data_list):
    y_min = -0.6
    y_max = 0.6
    for i in range(len(data_list)):
        theta1 = data_list[i][:,2]
        plt.plot(x, theta1, label= f"k={k_list[i]}")
        plt.vlines(x[horizon_entry[i]], color = colors[i], ymin=y_min, ymax=y_max, linestyle="dashdot")
    plt.title(r"$\Theta_1$")
    plt.xlabel("x")
    plt.xlim(-12,-0)
    plt.ylim(y_min,y_max)
    plt.vlines(x[g_tilde_max], ymin=y_min, ymax=y_max, linestyle="--", label=r"$x_{rec}$")
    #plt.vlines(x[Xe_half], ymin=-0.5, ymax=0.5, linestyle="--", label=r"x_{rec}")
    plt.axvspan(x[0],x[rad_matter_eq],color="peachpuff")
    plt.axvspan(x[rad_matter_eq], x[matter_de_eq], color="powderblue")
    plt.axvspan(x[matter_de_eq], x[-1], color="pink")
    #plt.legend(["k=0.1", "k=0.01", "k=0.001"])
    plt.legend(loc='lower right')
    plt.savefig("theta1_plot.pdf")
    plt.show()

plot_theta1(x, data_list)


def plot_theta2(x, data_list):
    y_min = -0.06
    y_max = 0.12    
    for i in range(len(data_list)):
        theta2 = data_list[i][:,3]
        plt.plot(x, theta2, color = colors[i], label= f"k={k_list[i]}")
        plt.vlines(x[horizon_entry[i]], color = colors[i], ymin=y_min, ymax=y_max, linestyle="dashdot")
    plt.title(r"$\Theta_2$")
    plt.xlabel("x")
    plt.xlim(-12,-0)
    plt.ylim(y_min,y_max)
    plt.vlines(x[g_tilde_max], ymin=y_min, ymax=y_max, linestyle="--", label=r"$x_{rec}$")
    #plt.vlines(x[Xe_half], ymin=-0.05, ymax=0.1, linestyle="--", label=r"x_{rec}")
    plt.axvspan(x[0],x[rad_matter_eq],color="peachpuff")
    plt.axvspan(x[rad_matter_eq], x[matter_de_eq], color="powderblue")
    plt.axvspan(x[matter_de_eq], x[-1], color="pink")
    #plt.legend(["k=0.1", "k=0.01", "k=0.001"])
    plt.legend()
    plt.savefig("theta2_plot.pdf")
    plt.show()

plot_theta2(x, data_list)



def plot_Phi(x, data_list):
    y_min = -0.2
    y_max = 0.8
    for i in range(len(data_list)):
        Phi = data_list[i][:,4]
        plt.plot(x, Phi, color = colors[i], label= f"k={k_list[i]}")
        plt.vlines(x[horizon_entry[i]], color = colors[i], ymin=y_min, ymax=y_max, linestyle="dashdot")
    plt.title(r"$\Phi$")
    plt.xlabel("x")
    plt.xlim(-12,-0)
    plt.ylim(y_min,y_max)
    plt.vlines(x[g_tilde_max], ymin=y_min, ymax=y_max, linestyle="--", label=r"$x_{rec}$")
    #plt.vlines(x[Xe_half], ymin=0, ymax=0.7, linestyle="--", label=r"x_{rec}")
    plt.axvspan(x[0],x[rad_matter_eq],color="peachpuff")
    plt.axvspan(x[rad_matter_eq], x[matter_de_eq], color="powderblue")
    plt.axvspan(x[matter_de_eq], x[-1], color="pink")
    #plt.legend(["k=0.1", "k=0.01", "k=0.001"])
    plt.legend(loc="lower right")
    plt.savefig("Phi_plot.pdf")
    plt.show()

plot_Phi(x, data_list)



def plot_Psi(x, data_list):
    y_min = -0.7
    y_max = 0.2
    for i in range(len(data_list)):
        Psi = data_list[i][:,5]
        plt.plot(x, Psi, color = colors[i], label= f"k={k_list[i]}")
        plt.vlines(x[horizon_entry[i]], color = colors[i], ymin=y_min, ymax=y_max, linestyle="dashdot")
    plt.title(r"$\Psi$")
    plt.xlabel("x")
    plt.xlim(-12,-0)
    plt.ylim(y_min,y_max)
    plt.vlines(x[g_tilde_max], ymin=y_min, ymax=y_max, linestyle="--", label=r"$x_{rec}$")
    #plt.vlines(x[Xe_half], ymin=-0.7, ymax=0, linestyle="--", label=r"x_{rec}")
    plt.axvspan(x[0],x[rad_matter_eq],color="peachpuff")
    plt.axvspan(x[rad_matter_eq], x[matter_de_eq], color="powderblue")
    plt.axvspan(x[matter_de_eq], x[-1], color="pink")
    #plt.legend(["k=0.1", "k=0.01", "k=0.001"])
    plt.legend(loc="lower right")
    plt.savefig("Psi_plot.pdf")
    plt.show()

plot_Psi(x, data_list)

def plot_vb_v_cdm(x, data_list):
    y_min = 1e-4
    y_max = 1.5e1
    for i in range(len(data_list)):
        vb = abs(data_list[i][:,6])
        vb_cdm = abs(data_list[i][:,7])
        plt.semilogy(x, vb,"--" , color = colors[i], label=r"$v_b$")
        plt.semilogy(x,vb_cdm, color = colors[i], label=r"$v_{CDM}$")
        plt.vlines(x[horizon_entry[i]], color = colors[i], ymin=y_min, ymax=y_max, linestyle="dashdot")
    plt.title(r"$|v_b|, |v_{CDM}$|")
    plt.xlabel("x")
    plt.xlim(-12,-0)
    plt.ylim(y_min,y_max)
    plt.vlines(x[g_tilde_max], ymin=y_min, ymax=y_max, linestyle="--", label=r"$x_{rec}$")
    #plt.vlines(x[Xe_half], ymin=1e-4, ymax=1.5e1, linestyle="--", label=r"x_{rec}")
    plt.axvspan(x[0],x[rad_matter_eq],color="peachpuff")
    plt.axvspan(x[rad_matter_eq], x[matter_de_eq], color="powderblue")
    plt.axvspan(x[matter_de_eq], x[-1], color="pink")
    #plt.legend(["k=0.1", "k=0.01", "k=0.001"])
    #plt.legend()
    plt.legend([r"$v_b$",r"$v_{cdm}$"])
    plt.savefig("v_plot.pdf")
    plt.show()

plot_vb_v_cdm(x, data_list)

def plot_deltas(x, data_list):
    y_min = 1e-4
    y_max = 1e5
    for i in range(len(data_list)):
        deltab = abs(data_list[i][:,8])
        delta_cdm = abs(data_list[i][:,9] )
        plt.semilogy(x, deltab,"--", color = colors[i],label=r"$\delta_b$")
        plt.semilogy(x,delta_cdm, color = colors[i],label=r"$\delta_{cdm}$")
        plt.vlines(x[horizon_entry[i]], color = colors[i], ymin=y_min, ymax=y_max, linestyle="dashdot")
    plt.title(r"$|\delta_b|, |\delta_{CDM}$|")
    plt.xlabel("x")
    plt.xlim(-12,-0)
    plt.ylim(y_min,y_max)
    plt.vlines(x[g_tilde_max], ymin=y_min, ymax=y_max, linestyle="--", label=r"$x_{rec}$")
    #plt.vlines(x[Xe_half], ymin=1e-4, ymax=1e5, linestyle="--", label=r"x_{rec}")
    plt.axvspan(x[0],x[rad_matter_eq],color="peachpuff")
    plt.axvspan(x[rad_matter_eq], x[matter_de_eq], color="powderblue")
    plt.axvspan(x[matter_de_eq], x[-1], color="pink")
    #plt.legend(["k=0.1", "k=0.01", "k=0.001"])
    plt.legend([r"$\delta_b$",r"$\delta_{cdm}$"])
    plt.savefig("delta_plot.pdf")
    plt.show()
        

plot_deltas(x, data_list)

