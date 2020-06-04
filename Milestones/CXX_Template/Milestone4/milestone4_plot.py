import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 14})
import os




cells = np.genfromtxt("../cells.txt")


n = int(len(cells[:,0]))
ells = cells[:,0]
C_ell = cells[:,0]

plt.semilogy(ells, C_ell)
plt.show()