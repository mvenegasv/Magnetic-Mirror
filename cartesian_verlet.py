import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
M = 1.0  # Normalized mass
Q = 1.0  # Normalized charge
q_m = Q / M  # Charge-to-mass ratio
Bz0 = 1.0  # Magnetic field at z = 0
dBz_dz = 1.0  # Magnetic field gradient

# Initial conditions
v0x = 0.2 #0.2
v0y = 0 #0.0
v0z = 0.3 #0.3
x0, y0, z0 = 0.1, 0.1,0.5  # Initial position (normalized)

# Magnetic field components
def Bz(z):
    """Axial magnetic field."""
    return Bz0 + z * dBz_dz

def Br(r):
    """Radial magnetic field derived from divergence-free condition."""
    return -0.5 * r * dBz_dz

# Cartesian coordinates 
def theta(x, y):
    """Angle in cylindrical coordinates."""
    return np.arctan2(y, x)

def Bx(x, y, z):
    """X-component of the magnetic field."""
    r = np.sqrt(x**2 + y**2)
    return Br(r) * np.cos(theta(x, y))

def By(x, y, z):
    """Y-component of the magnetic field."""
    r = np.sqrt(x**2 + y**2)
    return Br(r) * np.sin(theta(x, y))

# Verlet integration method
def verlet(state, dt, n_steps):
    x, y, z, vx, vy, vz = state
    positions = []
    velocities = []
    
    # Compute initial acceleration
    Bx_val = Bx(x, y, z)
    By_val = By(x, y, z)
    Bz_val = Bz(z)
    ax = q_m * (vy * Bz_val - vz * By_val)
    ay = q_m * (vz * Bx_val - vx * Bz_val)
    az = q_m * (vx * By_val - vy * Bx_val)
    
    for _ in range(n_steps):
        # Store positions and velocities
        positions.append([x, y, z])
        velocities.append([vx, vy, vz])

        # Update position
        x_new = x + vx * dt + 0.5 * ax * dt**2
        y_new = y + vy * dt + 0.5 * ay * dt**2
        z_new = z + vz * dt + 0.5 * az * dt**2

        # Compute magnetic field at the new position
        Bx_val_new = Bx(x_new, y_new, z_new)
        By_val_new = By(x_new, y_new, z_new)
        Bz_val_new = Bz(z_new)

        # Compute new acceleration
        ax_new = q_m * (vy * Bz_val_new - vz * By_val_new)
        ay_new = q_m * (vz * Bx_val_new - vx * Bz_val_new)
        az_new = q_m * (vx * By_val_new - vy * Bx_val_new)

        # Update velocity
        vx += 0.5 * (ax + ax_new) * dt
        vy += 0.5 * (ay + ay_new) * dt
        vz += 0.5 * (az + az_new) * dt

        # Update positions and accelerations for the next iteration
        x, y, z = x_new, y_new, z_new
        ax, ay, az = ax_new, ay_new, az_new

    return np.array(positions), np.array(velocities)

# Initial state: [x, y, z, vx, vy, vz]
state0 = [x0, y0, z0, v0x, v0y, v0z]

# Time step and number of steps
dt = 1e-4
n_steps = 350000#300000

# Solve using Verlet method
positions, velocities = verlet(state0, dt, n_steps)
x, y, z = positions[:, 0], positions[:, 1], positions[:, 2]
vx, vy, vz = velocities[:, 0], velocities[:, 1], velocities[:, 2]
time = np.linspace(0, dt * n_steps, n_steps)
radius = np.sqrt(x**2 + y**2)
v_perp2 = vx**2 + vy**2  # Perpendicular velocity squared
B_mag =np.sqrt(Bx(x, y, z)**2 + By(x, y, z)**2 + Bz(z)**2) # Magnetic field strength

# Check stability
B_max= max(B_mag)
omega_max= q_m * B_max
T_min= 2*np.pi/omega_max
print('smallest characteristic time=',T_min,'     recomended dt=',T_min*0.01)

# Plot the 3D trajectory
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, color='red', lw=2)
ax.set_title("Particle Trajectory")
ax.set_xlabel("x [m]" )
ax.set_ylabel("y [m] ")
ax.set_zlabel("z [m]")
plt.show()



#-------- Length scales of the system------------------
# Larmor radius
r_L = np.sqrt(v_perp2) / B_mag
 # Gradients in cylindrical coordinates
grad_Bz = dBz_dz  # Axial field gradient
grad_Br = - 0.5 * dBz_dz 
grad_Bx = grad_Br * np.cos(theta(x, y))
grad_By =grad_Br *  np.sin(theta(x, y))
grad_B = np.sqrt(grad_Bx**2 + grad_By**2 + grad_Bz**2)
 # Field gradient scale length
L_B = B_mag / grad_B
# Plot radial distance over time
plt.figure()
plt.plot(time, radius,color='black', lw=2)
plt.plot(time, r_L ,color='orange', lw=2,label='Larmor radius [m]')
plt.title("Radial Distance Over Time")
plt.xlabel("Time [s]" )
plt.ylabel("Radius [m] ")
plt.legend()
plt.show()

#Adiabatic Parameter
epsilon_values = r_L / L_B 
# Plot epsilon over time
plt.figure(figsize=(7, 3))
plt.plot(time, epsilon_values,color='black')
plt.title("Adiabatic Parameter ($\\epsilon$) Over Time")
plt.xlabel("Time [s]")
plt.ylabel("$\\epsilon$ ")
#plt.savefig('/Users/mv58/Documents/PhD/Plasma/Mirror25/epsilon_verlet.png', bbox_inches='tight')
plt.show()



# Plot height distance/ velocities over time
fig = plt.figure(figsize=(11, 4))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1])  
# Left subplot: 2D magnetic field with particle trajectory
ax1 = fig.add_subplot(gs[0])
ax1.plot(time, z,color='black', lw=2)
ax1.set_xlabel("Time [s]",size=15)
ax1.set_ylabel("Height [m]",size=15)
ax1.tick_params(axis='both', which='major', labelsize=15) 
ax2 = fig.add_subplot(gs[1])
ax2.plot(time, np.abs(vz),label='$v_\\parallel$')
ax2.plot(time, np.sqrt(vx**2 + vy**2),label='$v_\\perp$')
ax2.tick_params(axis='both', which='major', labelsize=15) 
ax2.set_xlabel("Time [s]",size=15)
ax2.set_ylabel("Velocities [m/s]",size=15)
plt.legend(fontsize=11,loc=3,frameon=False)
#plt.savefig('/Users/mv58/Documents/PhD/Plasma/Mirror25/height_velocities_verlet.png', dpi=300, bbox_inches='tight')
plt.show()



# magnetic moment
mu = v_perp2 / (2 * B_mag)
# Energy
v_total2 = vx**2 + vy**2 + vz**2  # Total velocity squared
energy = 0.5 * v_total2  # Kinetic energy (normalized)

# Plot total energy over time
fig, ax = plt.subplots()
ax.set_yscale('log')
ax.tick_params(axis='y',which='both', right=True, labelright=True)
ax.set_xlim(0,max(time))
ax.set_ylim(1e-2,1)
ax.plot(time, mu,color='red',label='magnetic moment [J/T]')
ax.plot(time, energy,color='blue',label='Total energy [J]' )
ax.set_title("Conserved quantities")
ax.set_xlabel("Time")
plt.legend()
plt.grid()
plt.show()




