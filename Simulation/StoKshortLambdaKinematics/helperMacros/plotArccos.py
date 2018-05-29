import numpy as np
import matplotlib.pyplot as plt

def f(x):
        return np.arccos(np.sqrt(x-1)/np.sqrt(3))
def g(x):
        return np.arccos(-np.sqrt(x-1)/np.sqrt(3))


x1 = np.arange(0.0, 5.0, 0.1)

plt.plot(x1, f(x1), 'bo')
plt.plot(x1, g(x1), 'ro')
plt.plot(x1, 2*g(x1)-3.1415, 'go')

plt.show()
