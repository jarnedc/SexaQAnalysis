import matplotlib.pyplot as plt
import numpy as np

M_S=1.8
M_S2=2
q=1.15
T=0.08

x = np.arange(0.0, 10.0, 0.01)

def TSalis(x,M_S,q,T):
    y = 10*x*np.sqrt(x*x+M_S*M_S)*np.power(1+(q-1)*np.sqrt(x*x+M_S*M_S)/T,-1.15/(q-1))
    return y

y1 = TSalis(x,1.8,q,T)
y2 = TSalis(x,1.8,q,T+0.01)
y3 = TSalis(x,1.8,q,T-0.01)

plt.plot(x,y1,'b-',x,y2,'g-',x,y3,'r-')
plt.show()
