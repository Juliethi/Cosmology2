import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 14})
import os
import seaborn as sn
sn.set_style("darkgrid")


Mpc = 3.08567758e22
h = 0.7
T_cmb = 2.725

cells_full = np.genfromtxt("../cells_full.txt")
cells_nolastterm = np.genfromtxt("../cells_nolastterm.txt")

cells_sw = np.genfromtxt("../cells_SW.txt")
cells_isw = np.genfromtxt("../cells_ISW.txt")
cells_doppler = np.genfromtxt("../cells_Doppler.txt")
cells_lastterm = np.genfromtxt("../cells_lastterm.txt")

ells = cells_full[:,0]


"""DIFFERENT TEMPERATURE POWER SPECTRUM PLOTS"""

def plot_power_spectrum_full():
    cells = cells_full[:,1]#/(1e6*T_cmb)**2
    plt.loglog(ells, cells)
    plt.title(r"$C_{\ell}$ full")
    plt.xlabel(r"$\ell$")
    plt.ylabel(r"$\frac{\ell(\ell+1)}{2\pi}(10^6 T_{\mathrm{CMB}})^2$")
    plt.tight_layout()
    plt.savefig("cells_full.pdf")
    plt.show()
    return None

def plot_power_spectrum_components():
    plt.loglog(ells, cells_full[:,1], c="k")
    plt.loglog(ells, cells_sw[:,1])
    plt.loglog(ells, cells_isw[:,1])
    plt.loglog(ells, cells_doppler[:,1])
    plt.loglog(ells, cells_lastterm[:,1])
    #plt.ylim(1e-1,1e4)
    plt.title(r"$C_{\ell}$ for each source function term")
    plt.xlabel(r"$\ell$")
    plt.ylabel(r"$\frac{\ell(\ell+1)}{2\pi}(10^6 T_{\mathrm{CMB}})^2 C_{\ell}$")
    plt.tight_layout()
    plt.legend(["full","SW", "ISW", "Doppler", "Last term"])
    plt.savefig("cells_component.pdf")
    plt.show()
    return None

def plot_power_spectrum_compare():
    plt.loglog(ells, cells_full[:,1], label="All terms")
    plt.loglog(ells, cells_nolastterm[:,1], label="No last term")
    plt.title(r"$C_{\ell}$ comparison")
    plt.xlabel(r"$\ell$")
    plt.ylabel(r"$\frac{\ell(\ell+1)}{2\pi}(10^6 T_{\mathrm{CMB}})^2 C_{\ell}$")
    plt.legend()
    plt.tight_layout()
    plt.savefig("cells_compare.pdf")
    plt.show()
    return None

def plot_power_spectrum_nolastterm():
    plt.loglog(ells, cells_nolastterm[:,1])
    plt.title(r"$C_{\ell}$ without last source function term")
    plt.xlabel(r"$\ell$")
    plt.ylabel(r"$\frac{\ell(\ell+1)}{2\pi}(10^6 T_{\mathrm{CMB}})^2 C_{\ell}$")
    plt.tight_layout()
    plt.savefig("cells_nolastterm.pdf")
    plt.show()
    return None

#plot_power_spectrum_full()
#plot_power_spectrum_components()
#plot_power_spectrum_compare()
#plot_power_spectrum_nolastterm()

"""PLOT THETAS"""
ells_list = [
        2,    3,    4,    5,    6,    7,    8,    10,   12,   15,   
        20,   25,   30,   40,   50,   60,   70,   80,   90,   100,  
        120,  140,  160,  180,  200,  225,  250,  275,  300,  350,  
        400,  450,  500,  550,  600,  650,  700,  750,  800,  850,  
        900,  950,  1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 
        1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 
        1900, 1950, 2000]

iell_list = [0,10,20,60]

theta_list = []

for i in range(len(iell_list)):
    theta_list.append(np.genfromtxt(f"../thetas_{iell_list[i]}.txt"))
k_theta = theta_list[0][:,0]*Mpc

def plot_transfer_function():
    for i in range(len(theta_list)):
        plt.loglog(k_theta, theta_list[i][:,1], label = f"$\ell$={ells_list[iell_list[i]]}")
    plt.yscale("symlog", linthreshy=1e-5)
    plt.title(r"Transfer function $\Theta_{\ell}(k)$")
    plt.xlabel(r"k [1/Mpc]")
    plt.ylabel(r"$\Theta_{\ell}(k)$")
    plt.legend()
    plt.tight_layout()
    plt.savefig("transerfunction.pdf")
    plt.show()
    return None

def plot_theta_integrand():
    for i in range(len(theta_list)):
        integrand = theta_list[i][:,1]**2/k_theta
        plt.loglog(k_theta, integrand, label = f"$\ell$={ells_list[iell_list[i]]}")
    #plt.yscale("symlog", linthreshy=1e-5)
    plt.title(r"Spectrum integrand $\Theta_{\ell}(k)^2/k$")
    plt.xlabel(r"k [1/Mpc]")
    plt.ylabel(r"$\Theta_{\ell}(k)^2/k$")
    plt.ylim([1e-12,1e2])
    plt.legend()
    plt.tight_layout()
    plt.savefig("thetaintegrand.pdf")
    plt.show()   


#plot_transfer_function()
plot_theta_integrand()

"""PLOT MATTER POWER SPECTRUM"""

matter = np.genfromtxt("../powerspectrum.txt")
k_matter = matter[1:,0]*Mpc/h
pk = matter[1:,1]*(h/Mpc)**3

k_peak = matter[0,0]*Mpc/h
print(k_peak)

def plot_matter_powerspectrum():
    plt.loglog(k_matter, pk, label="P(k)")
    plt.title("Matter power spectrum")
    plt.vlines(k_peak, pk.min(), pk.max(), label=r"$k_{peak}$", colors = "r")
    plt.xlabel("k [h/Mpc]")
    plt.ylabel(r"P(k) [Mpc/h]$^3$")
    plt.ylim(pk.min(), 0.5*1e5)
    plt.legend()
    plt.tight_layout()
    plt.savefig("matterpowerspectrum.pdf")
    plt.show()

plot_matter_powerspectrum()

