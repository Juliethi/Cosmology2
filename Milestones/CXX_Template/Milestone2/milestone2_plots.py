import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 14})

values = np.genfromtxt("../recombination.txt")
#x, Xe(X), ne(X), tau,tau', tau'', g,g',g''

n = len(values[:,0])
x = values[:,0]
z = 1/np.exp(x) - 1

Xe = values[:,1]

plt.semilogy(x,Xe)
plt.show()