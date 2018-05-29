import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D

M_S=2
M_S2=2
q=1.15
T=0.08


def PZTSalisNorm(pT,M_S,q,T,eta):
    y = (pT*np.power(1+(q-1)*np.sqrt(pT*pT+M_S*M_S)/T,-q/(q-1)))/integrate.quad(lambda x : x*np.power(1+(q-1)*np.sqrt(x*x+M_S*M_S)/T,-q/(q-1)),0,10000)[0] / np.tan(2*np.arctan(np.exp(-eta)))
    return y

r_eta = np.linspace (-2.5, 2.5, 100)
r_pt = np.linspace(0,10, 100)

R_eta, R_pt = np.meshgrid(r_eta, r_pt)
r_z = PZTSalisNorm(R_pt,M_S,q,T,R_eta)

fig = plt.figure()
ax = plt.axes(projection='3d')
#ax.contour3D(np.abs(R_eta), R_pt, r_z, 50, cmap='binary')
ax.plot_wireframe(np.abs(R_eta), R_pt, r_z, rstride=10, cstride=10)
ax.plot_surface(np.abs(R_eta), R_pt, r_z, color = 'b')
ax.set_xlabel('abs(eta)')
ax.set_ylabel('pt')
ax.set_zlabel('Pz');

plt.title("S pz TSalis distribution ")
fig = plt.gcf()
plt.show()
pp = PdfPages('TSalisPz.pdf')
pp.savefig(fig)
pp.close()

