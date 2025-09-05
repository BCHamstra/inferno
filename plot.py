import numpy as np
import matplotlib.pyplot as plt
v = np.array([5e-2, 1e-1, 2e-1, 3e-1, 4e-1, 6e-1, 8e-1])
c = np.array([0.75, 0.55, 0.4,  0.3,  0.25, 0.2,  0.15])
h = np.array([2.1,  2.3,  2.45, 2.55, 2.55, 2.65, 2.65])

vinv = 1/v

plt.plot(vinv, 3-h, label="H")
plt.plot(vinv, c, label="C")
plt.xlabel("1/v [s/cm]")
plt.ylabel("$\Delta y$ [cm]")
plt.legend()
plt.show()
