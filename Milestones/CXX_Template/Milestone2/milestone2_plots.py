import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
sns.set_style("darkgrid")
matplotlib.rcParams.update({'font.size': 14})

values = np.genfromtxt("../recombination.txt")
#x, Xe(X), Xe_Saha_only, ne(X), tau,tau', tau'', g,g',g''

n = len(values[:,0])
x = values[:,0]
z = 1/np.exp(x) - 1

Xe = values[:,1]
Xe_Saha = values[:,2]

tau = values[:,4]
dtau = values[:,5]
ddtau = values[:,6]

g = values[:,7]
dg = values[:,8]
ddg = values[:,9]

def find_xstar_zstar(x, z, tau):
    tau_diff = abs(tau-1.)
    t1 = np.argmin(tau_diff) #Finding where abs(tau-1) is closes to 0, that is where tau is approximately 1
    x_star = x[t1]
    z_star = z[t1]
    print("x when tau(x)=1 is = ", x_star)
    print("z when tau(z)=1 is =", z_star)
    return x_star

x_star = find_xstar_zstar(x, z, tau)

def find_xrec_zrec(x, z, Xe, Xe_Saha):
    Xe_diff = abs(Xe-0.5)
    Xe_Saha_diff = abs(Xe_Saha-0.5)
    Xe05 = np.argmin(Xe_diff)
    Xe_Saha05 = np.argmin(Xe_Saha_diff)
    
    x_rec = x[Xe05]; z_rec = z[Xe05]
    x_rec_saha = x[Xe_Saha05]; z_rec_saha = z[Xe_Saha05]

    print("x and z values when Xe = 0.5", x_rec, z_rec)
    print("x and z values when Xe from Saha only = 0.5", x_rec_saha, z_rec_saha)
    return x_rec, x_rec_saha

x_rec, x_rec_saha = find_xrec_zrec(x, z, Xe, Xe_Saha)


def plot_xe(x, Xe, Xe_Saha):
    fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(14,6))
    ax1.semilogy(x,Xe, label=r"$X_e$")
    ax1.semilogy(x, Xe_Saha, label=r"$X_e$ (Saha)")
    ax1.set_xlim(-12,-0)
    ax1.set_ylim(1e-4,0.3*1e1)
    ax1.set_xlabel("x")
    ax1.set_ylabel(r"$X_e(x)$")
    ax1.legend()
    
    ax2.semilogy(x,Xe, label=r"$X_e$")
    ax2.semilogy(x, Xe_Saha, label=r"$X_e$ (Saha)")
    ax2.set_xlim(-7.6,-6.8)
    ax2.set_ylim(1e-1,0.3*1e1)
    ax2.vlines(x_rec, ymin=1e-4, ymax=3*1e1, colors="g", linestyle="-.",label=r"$x_{rec}$")
    ax2.vlines(x_rec_saha, ymin=1e-4, ymax=3*1e1, colors="y", linestyle="-." ,label=r"$x_{rec}$ (Saha)")
    ax2.set_xlabel("x")
    ax2.legend()
    
    fig.tight_layout()
    fig.savefig("xe.pdf")
    plt.show()

def plot_tau(x, tau, dtau, ddtau):
    plt.semilogy(x,tau,label=r"$\tau$")
    plt.semilogy(x,-dtau, label=r"$-\tau$'")
    plt.semilogy(x,ddtau, label=r"$\tau$''")
    plt.title(r"$\tau(x)$")
    plt.ylim(1e-8,1e6)
    plt.xlim(-12,-0)
    plt.vlines(x_star, 1e-8, 1e6, color="mediumpurple", linestyle="-.", label=r"$x_{*}$")
    plt.legend()
    plt.savefig("optical_depth.pdf")
    plt.show()

def plot_g(x, g, dg, ddg):
    plt.plot(x,g, label=r"$\tilde{g}$")
    plt.plot(x,dg, label=r"$\tilde{g}$'")
    plt.plot(x,ddg, label=r"$\tilde{g}$''")
    plt.title(r"Unscaled visibility function $\tilde{g}$")
    plt.ylabel(r"$\tilde{g}(x)$")
    plt.xlabel("x")
    plt.ylim(-700,700)
    plt.xlim(-12,-0)
    plt.legend()
    plt.savefig("g.pdf")
    plt.show()

def plot_g_scaled(x, g, dg, ddg):
    plt.plot(x,g, label=r"$\tilde{g}$")
    plt.plot(x,dg/10., label=r"$\tilde{g}'/10$")
    plt.plot(x,ddg/100., label=r"$\tilde{g}''/100$")
    plt.title(r"Scaled visibility function $\tilde{g}$")
    plt.xlim(-8,-6)
    plt.ylim(-9,9)
    plt.ylabel(r"$\tilde{g}(x)$")
    plt.xlabel("x")
    plt.vlines(x_star, -9, 9, color="mediumpurple", linestyle="-.", label=r"$x_{*}$")
    plt.legend()
    plt.savefig("g_scaled.pdf")
    plt.show()


plot_xe(x, Xe, Xe_Saha)
plot_tau(x, tau, dtau, ddtau)
plot_g(x, g, dg, ddg)
plot_g_scaled(x, g, dg, ddg)