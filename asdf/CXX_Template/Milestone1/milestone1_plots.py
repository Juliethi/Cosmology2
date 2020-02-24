import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import seaborn as sns
matplotlib.rcParams.update({'font.size': 14})
#sns.set()

values = np.genfromtxt("../cosmology.txt")
#x, eta(x), Hp(x), dHp_dx, omegaB, omegaCDM, omegaLambda, OmegaR, OmegaNu, OmegaK"

n = len(values[:,0])
Mpc = 3.08567758e22 #Mpc in meters
x = values[:,0]
z = 1/np.exp(x) - 1

eta = values[:,1]/Mpc/1000. #Converting eta from m to Mpc

Hp = values[:,2]*Mpc/1000. #Converting Hp to km/s/Mpc
H = Hp/np.exp(x)

omegaB = values[:,4]
omegaCDM = values[:,5]
omegaLambda = values[:,6]
omegaR = values[:,7]

omegaSum = omegaB + omegaCDM + omegaLambda + omegaR
omegaMatter = omegaCDM + omegaB

def find_eras():
    rad_matter_dif = abs(omegaMatter - omegaR)
    rad_matter_dif = rad_matter_dif[:int(n/2)]
    rad_matter_eq = np.argmin(rad_matter_dif)
    print("x value rad matter eq:", x[rad_matter_eq])
    print("z value of rad matter eq:", z[rad_matter_eq])
    
    matter_de_dif = abs(omegaMatter - omegaLambda)
    matter_de_dif = matter_de_dif[int(n/2):]
    matter_de_eq = np.argmin(matter_de_dif) + int(n/2)
    print("x value matter de eq:", x[matter_de_eq])
    print("z value of matter de eq:", z[matter_de_eq])

    return rad_matter_eq, matter_de_eq

rad_matter_eq, matter_de_eq = find_eras()
print(rad_matter_eq, matter_de_eq)

def plot_omegas():
    plt.plot(x, omegaLambda, label=r"$\Omega_{\Lambda}$")
    plt.plot(x, omegaR, label=r"$\Omega_R$")
    plt.plot(x, (omegaB + omegaCDM), label=r"$\Omega_B + \Omega_{CDM}$")
    plt.plot(x, omegaSum, "-.", label=r"$\sum \Omega_i$")
    plt.plot(x,omegaB, "--", label=r"$\Omega_B$")
    plt.plot(x,omegaCDM, "--",label=r"$\Omega_{CDM}$")
    plt.axvspan(x[0],x[rad_matter_eq],color="peachpuff")
    plt.axvspan(x[rad_matter_eq], x[matter_de_eq], color="powderblue")
    plt.axvspan(x[matter_de_eq], x[-1], color="pink")
    plt.legend()
    plt.title(r"Time evolution of $\Omega_i$")
    plt.xlabel("x")
    plt.ylabel(r"$\Omega$")
    plt.xlim(x[0],x[-1])
    plt.savefig("omega_plot.pdf")
    plt.show()

#plot_omegas()

def plot_H_eta():
    matplotlib.rcParams.update({'font.size': 20})
    fig, ax = plt.subplots(2,2,figsize=(14,10))
    
    ax[0,0].semilogy(x, H)
    ax[0,0].set_title("H(x) [km/s/Mpc]")
    ax[0,0].set_xlabel("x")
    ax[0,0].set_ylabel("H")
    ax[0,0].axvspan(x[0],x[rad_matter_eq],color="peachpuff")
    ax[0,0].axvspan(x[rad_matter_eq], x[matter_de_eq], color="powderblue")
    ax[0,0].axvspan(x[matter_de_eq], x[-1], color="pink")
    ax[0,0].set_xlim(x[0],x[-1])

    ax[0,1].loglog(z, H)
    ax[0,1].set_title("H(z) [km/s/Mpc]")
    ax[0,1].set_xlabel("z")
    ax[0,1].set_ylabel("H")
    #ax[0,1].set_xlim(z[0],z[-1])
    ax[0,1].invert_xaxis()
    ax[0,1].axvspan(z[0],z[rad_matter_eq],color="peachpuff")
    ax[0,1].axvspan(z[rad_matter_eq], z[matter_de_eq], color="powderblue")
    ax[0,1].axvspan(z[matter_de_eq], z[-1], color="pink")
    ax[0,1].set_xlim(z[0],z[-1])

    ax[1,0].semilogy(x, Hp)
    ax[1,0].set_title(r"$\mathcal{H}(x)$ [km/s/Mpc]")
    ax[1,0].set_xlabel("x")
    ax[1,0].set_ylabel(r"$\mathcal{H}$")
    ax[1,0].axvspan(x[0],x[rad_matter_eq],color="peachpuff")
    ax[1,0].axvspan(x[rad_matter_eq], x[matter_de_eq], color="powderblue")
    ax[1,0].axvspan(x[matter_de_eq], x[-1], color="pink")
    ax[1,0].set_xlim(x[0],x[-1])


    ax[1,1].semilogy(x, eta)
    ax[1,1].set_title(r"$\eta(x)$ [Gpc]")
    ax[1,1].set_xlabel("x")
    ax[1,1].set_ylabel(r"$\eta$")
    ax[1,1].axvspan(x[0],x[rad_matter_eq],color="peachpuff")
    ax[1,1].axvspan(x[rad_matter_eq], x[matter_de_eq], color="powderblue")
    ax[1,1].axvspan(x[matter_de_eq], x[-1], color="pink")
    ax[1,1].set_xlim(x[0],x[-1])

    fig.tight_layout()
    fig.savefig("H_eta_plots.pdf")
    plt.show()
    

plot_H_eta()
