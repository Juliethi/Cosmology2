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

tau = values[:,3]
dtau = values[:,4]
ddtau = values[:,5]

g = values[:,6]
dg = values[:,7]
ddg = values[:,8]



plt.semilogy(x,Xe)
plt.title("Xe")
plt.show()


plt.semilogy(x,tau,label="tau")
plt.semilogy(x,-dtau, label="dtau")
plt.semilogy(x,ddtau, label="ddtau")
plt.title("tau")
plt.ylim(1e-8,1e6)
plt.legend()
plt.show()


plt.plot(x,g, label="g")
plt.plot(x,dg, label="dg")
plt.plot(x,ddg, label="ddg")
plt.title("g")
plt.legend()
plt.show()