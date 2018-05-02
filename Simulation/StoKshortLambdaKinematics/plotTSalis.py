import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from matplotlib.backends.backend_pdf import PdfPages

M_S=2
M_S2=2
q=1.15
T=0.08

x1 = np.arange(0.0, 10.0, 0.01)

def TSalisNorm(x1,M_S,q,T):
    y = (x1*np.power(1+(q-1)*np.sqrt(x1*x1+M_S*M_S)/T,-q/(q-1)))/integrate.quad(lambda x : x*np.power(1+(q-1)*np.sqrt(x*x+M_S*M_S)/T,-q/(q-1)),0,10000)[0]
    #y = 10*x*np.sqrt(x*x+M_S*M_S)*np.power(1+(q-1)*np.sqrt(x*x+M_S*M_S)/T,-q/(q-1))
    return y


y1 = TSalisNorm(x1,M_S,q,     T)
y2 = TSalisNorm(x1,M_S,q+0.05,T)
y3 = TSalisNorm(x1,M_S,q+0.10,T)


plt.plot(x1,y1,'b-', label = 'q = 1.15, M_S = 2GeV, T = 0.08GeV')
plt.plot(x1,y2,'g-', label = 'q = 1.20, M_S = 2GeV, T = 0.08GeV')
plt.plot(x1,y3,'r-', label = 'q = 1.25, M_S = 2GeV, T = 0.08GeV')
plt.legend()
plt.xlabel("S pt")
plt.title("S pt TSalis distribution ")
fig = plt.gcf()
plt.show()
pp = PdfPages('TSalis.pdf')
pp.savefig(fig)
pp.close()

