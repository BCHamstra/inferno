import numpy as np
import matplotlib.pyplot as plt
# bit more compact than I would like to fit within max lines
# model parameters
nx  = 600  # number of points along length
ny  = 60   # number of points along height
L   = 30   # length of channel [cm]
h   = 3    # height of channel [cm]
p_s = 2    # start positon of temperature plates [cm]
p_e = 4    # end position of temperature plates [cm]
s   = 14   # sensor position [cm]
T0  = 20   # starting temperature [degC]
Th  = 40   # temperature hot plate [degC]
Tc  = 0    # temperature cold plate [degC]
v   = 1e-1 # fluid velocity [cm/s]
thr = 1e-6 # threshold between old and new values to stop simulation

poiseuille = True # whether to use Pouseuille or constant flow

# liquid constants
rho = 1e-3     # density [kg/cm3]
cp  = 4184     # heat capacity [J/kg K]
la  = 6.089e-3 # thermal conductivity [W/cm K]
mu  = 1e-3     # dynamic viscosity [Pa s]

# setup
dx    = L/nx
dy    = h/ny
T     = T0*np.ones((nx,ny))
T_new = np.zeros((nx,ny))
dif   = 1
runs  = 0
first  = True
second = True
for i in range(nx):
    if first and i*dx > p_s:
        i_s = i - 1
        first = False
    if not first and second and i*dx > p_e:
        i_e = i # index itself as it is not included in slices
        second = False
    if not second and i*dx > s:
        i_sens = i - 1
        break
x = np.linspace(0, L, num=nx, endpoint=False)
y = np.linspace(0, h, num=ny, endpoint=False)
if poiseuille:
    vx = 6*v*y/h*(1 - y/h)
else: # constant flow
    vx = v

# main loop
while dif > thr:
    # get temperature arrays for adjacent volumes
    T_yp = np.roll(T, (0,-1), axis=(0,1))
    T_ym = np.roll(T, (0,1), axis=(0,1))
    T_xm = np.roll(T, (1,0), axis=(0,1))
    T_xp = np.roll(T, (-1,0), axis=(0,1))
    # apply boundary conditions
    T_yp[:,-1]    = T[:,-1] # insulated wall: same temperature as fluid
    T_ym[:,0]   = T[:,0] 
    T_xp[-1,:] = T0 # fluid flow: far away has starting temperature
    T_xm[0,:]   = T0
    # boundary conditions for temperature plates
    T_yp[i_s:i_e,-1]  = Th
    T_ym[i_s:i_e,0] = Tc
    # determine new temperature
    T_new = (vx*rho*cp*T_xm*dy + la*((T_xp + T_xm)*dy/dx + (T_yp + T_ym)*dx/dy)) / (vx*rho*cp*dy + 2*la*(dx*dx+dy*dy)/(dx*dy))
    # determine difference
    dif = np.sum(np.abs(T_new - T))/(nx*ny)
    # update temperature
    T = np.copy(T_new)
    runs += 1

print("Runs: " + str(runs))
first = True
for i in range(ny):
    if first and T[i_sens,i] > 18:
        print("18 degrees at " + str((i-1) * dy) + " cm")
        first = False
    if not first and T[i_sens,i] > 22:
        print("22 degrees at " + str((i) * dy) + " cm")
        break
    
# plot result
pcm = plt.pcolormesh(np.repeat(x[:,np.newaxis], ny, axis=1), np.repeat(y[np.newaxis, :], nx, axis=0), T, cmap="coolwarm", shading="nearest")
cbar = plt.colorbar(pcm)
cbar.set_label("T [$\degree$C]")
plt.axvline(x=i_s*dx, color = "k", linestyle=":", label="Plates")
plt.axvline(x=i_e*dx, color = "k", linestyle=":")
plt.axvline(x=i_sens*dx, color = "g", linestyle=":", label='Sensor')
plt.legend()
plt.xlabel("x [cm]")
plt.ylabel("y [cm]")
plt.figure()
plt.plot(y, T[i_sens,:])
plt.xlabel("y [cm]")
plt.ylabel("T [$\degree$C]")
plt.show()
