import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 14})
import os




cells = np.genfromtxt("../cells.txt")
integrand = np.genfromtxt("../integrand.txt")

x = integrand[:,0]
pls = integrand[:,1]/1e-3

plt.plot(x, pls)
plt.xlim(-8,0)
plt.show()


n = int(len(cells[:,0]))
ells = cells[:,0]
C_ell = cells[:,1]

plt.loglog(ells, C_ell)
plt.show()

cells_sw = np.genfromtxt("../cells_SW.txt")
cells_isw = np.genfromtxt("../cells_ISW.txt")
cells_doppler = np.genfromtxt("../cells_Doppler.txt")
cells_lastterm = np.genfromtxt("../cells_lastterm.txt")

plt.loglog(ells, C_ell, c="k")
plt.loglog(ells, cells_sw[:,1])
plt.loglog(ells, cells_isw[:,1])
plt.loglog(ells, cells_doppler[:,1])
plt.loglog(ells, cells_lastterm[:,1])
plt.ylim(1e-1,1e4)
plt.legend(["full","sw", "isw", "Doppler", "D"])
#plt.show()