# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 10:19:51 2025

@author: mv58
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter
import matplotlib as mpl

# Set the correct path to FFmpeg
mpl.rcParams['animation.ffmpeg_path'] = '/opt/homebrew/bin/ffmpeg'


# Constants
M = 1.0  # Normalized mass
Q = 1.0  # Normalized charge
q_m = Q / M  # Charge-to-mass ratio
Bz0 = 1.0  # Magnetic field at z = 0
dBz_dz = 1.0  # Magnetic field gradient 

# Magnetic field components
# Cylindrical coordinates 
def Bz(z):
    """Axial magnetic field."""
    return Bz0 + z * dBz_dz

def Br(r):
    """Radial magnetic field derived from divergence-free condition."""
    return -0.5 * r * dBz_dz

def dBr_dr(r):
    """Radial gradient."""
    return -0.5 * dBz_dz


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

# Integrator
def leapfrog(state, dt, n_steps):
    x, y, z, vx, vy, vz = state
    positions = []
    velocities = []  # To store velocities

    for _ in range(n_steps):
        # Magnetic field at current position
        Bx_val = Bx(x, y, z)
        By_val = By(x, y, z)
        Bz_val = Bz(z)

        # Half-step velocity update
        vx_half = vx + 0.5 * dt * q_m * (vy * Bz_val - vz * By_val)
        vy_half = vy + 0.5 * dt * q_m * (vz * Bx_val - vx * Bz_val)
        vz_half = vz + 0.5 * dt * q_m * (vx * By_val - vy * Bx_val)

        # Full-step position update
        x += dt * vx_half
        y += dt * vy_half
        z += dt * vz_half

        # updating Magnetic field at the new postions
        Bx_val = Bx(x, y, z)
        By_val = By(x, y, z)
        Bz_val = Bz(z)

        # Full-step velocity update
        vx = vx_half + 0.5 * dt * q_m * (vy_half * Bz_val - vz_half * By_val)
        vy = vy_half + 0.5 * dt * q_m * (vz_half * Bx_val - vx_half * Bz_val)
        vz = vz_half + 0.5 * dt * q_m * (vx_half * By_val - vy_half * Bx_val)

        # Store positions and velocities
        positions.append([x, y, z])
        velocities.append([vx, vy, vz])

    return np.array(positions), np.array(velocities)

# Initial conditions
v0x = 0.2#0.1 #0.2
v0y = 0.
v0z = 0.3 #0.2 #0.4 #0.7
x0, y0, z0 = 0.1, 0.1, 0.5  # Initial position 

# Initial state: [x, y, z, vx, vy, vz]
state0 = [x0, y0, z0, v0x, v0y, v0z]

# Time step and number of steps
dt = 1e-4#2e-4 
n_steps = 350000 #360000 

# Solve using leapfrog method
positions, velocities = leapfrog(state0, dt, n_steps)

# Extract position and velocity components
x, y, z = positions[:, 0], positions[:, 1], positions[:, 2]
vx, vy, vz = velocities[:, 0], velocities[:, 1], velocities[:, 2]
time = np.linspace(0, dt * n_steps, n_steps)
radius = np.sqrt(x**2 + y**2)
v_perp2 = vx**2 + vy**2  # Perpendicular velocity squared
B_mag =np.sqrt(Bx(x, y, z)**2 + By(x, y, z)**2 + Bz(z)**2) # Magnetic field strength

    
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
plt.figure()
plt.plot(time, epsilon_values)
plt.title("Adiabatic Parameter ($\\epsilon$) Over Time")
plt.xlabel("Time ")
plt.ylabel("$\\epsilon$ ")
plt.grid()
plt.show()


# Plot height distance/ velocities over time
fig = plt.figure(figsize=(14, 6))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1])  
# Left subplot: 2D magnetic field with particle trajectory
ax1 = fig.add_subplot(gs[0])
ax1.plot(time, z,color='black', lw=2)
ax1.set_xlabel("Time [s]",size=15)
ax1.set_ylabel("Height [m]",size=15)
ax1.tick_params(axis='both', which='major', labelsize=15) 
ax2 = fig.add_subplot(gs[1])
ax2.plot(time, np.abs(vz),label='Parallel velocity')
ax2.plot(time, np.sqrt(vx**2 + vy**2),label='Perpendicular velocity')
ax2.tick_params(axis='both', which='major', labelsize=15) 
ax2.set_xlabel("Time [s]",size=15)
ax2.set_ylabel("Velocities [m/s]",size=15)
plt.legend(fontsize=13,loc=2)
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



#Adiabatic Parameter
epsilon_values = r_L / L_B 

# Plot epsilon over time
plt.figure()
plt.plot(time, epsilon_values)
plt.title("Adiabatic Parameter ($\\epsilon$) Over Time")
plt.xlabel("Time ")
plt.ylabel("$\\epsilon$ ")
plt.grid()
plt.show()



#%%  Animation 1
# Create a grid in Cartesian coordinates
x_grid = np.linspace(-0.03, 0.3, 100)  # x range (m)
z_grid = np.linspace(0, 3.5, 100)  # z range (m)
X, Z = np.meshgrid(x_grid, z_grid)  # 2D slice in the x-z plane
R = np.sqrt(X**2)  # Radial distance in the x-z plane

# Compute magnetic field components and magnitude
Bz_vals = Bz(Z)
Br_vals = Br(R)
Bx_vals = np.where(R != 0, Br_vals * X / R, 0)  
B_mag_vals = np.sqrt(Bx_vals**2 + Bz_vals**2)

# Create a combined figure with two subplots
fig = plt.figure(figsize=(12, 6))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1])  

# Left subplot: 2D magnetic field with particle trajectory
ax1 = fig.add_subplot(gs[0])
c = ax1.contourf(X, Z, np.log10(B_mag_vals), levels=50, cmap='plasma')
cb = fig.colorbar(c, ax=ax1, label='|B| [T]')
ax1.set_xlabel('x [m]')
ax1.set_ylabel('z [m]')
ax1.tick_params(axis='both', which='major', labelsize=6) 
ax1.set_title('Magnetic Field Strength |B| (x-z Plane)')
particle_2d, = ax1.plot([], [], 'ro')

# Right subplot: 3D particle trajectory
ax2 = fig.add_subplot(gs[1], projection='3d')
ax2.set_xlabel('x [m]')
ax2.set_ylabel('y [m]')
ax2.set_zlabel('z [m]')
ax2.set_title('3D Trajectory of the Particle')
ax2.tick_params(axis='both', which='major', labelsize=6)
line_3d, = ax2.plot([], [], [], lw=1, color='red')
point_3d, = ax2.plot([], [], [], 'ro')
time_text = ax2.text2D(0.05, 0.95, '', transform=ax2.transAxes)
ax2.set_xlim([np.min(x), np.max(x)])
ax2.set_ylim([np.min(y), np.max(y)])
ax2.set_zlim([np.min(z), np.max(z)])

# Downsample the trajectory to reduce the number of frames
step = 300 # Keep every 300 points
x1 = x[::step]
y1=y[::step]
z1 = z[::step]
time1=time[::step]

# Initialization function for the animation
def init():
    particle_2d.set_data([], [])
    line_3d.set_data([], [])
    line_3d.set_3d_properties([])
    point_3d.set_data([], [])
    point_3d.set_3d_properties([])
    time_text.set_text('')
    return particle_2d, line_3d, point_3d, time_text

# Update function for the animation
def update(frame):
    # Update 2D plot
    particle_2d.set_data([x1[frame]], [z1[frame]])
    
    # Update 3D plot
    line_3d.set_data(x1[:frame], y1[:frame])
    line_3d.set_3d_properties(z1[:frame])
    point_3d.set_data([x1[frame]], [y1[frame]])
    point_3d.set_3d_properties([z1[frame]])
    time_text.set_text(f'Time: {time1[frame]:.1f} s')
    
    return particle_2d, line_3d, point_3d, time_text

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=len(x1), init_func=init, blit=False, interval=20)

# Save the combined animation to a file
writer = FFMpegWriter(fps=50, metadata=dict(artist='Me'), bitrate=1800)
ani.save('/Users/mv58/Documents/PhD/Plasma/Mirror25/3dtrajectory_2.mp4', writer=writer)

plt.show()



#%%  Animation2:
fig = plt.figure()

# Conserved quantities
ax1 = fig.add_subplot()
ax1.set_xlabel('Time [s]')
ax1.set_yscale('log')
ax1.set_xlim(0, max(time1))
ax1.set_ylim(1e-2, 1)
ax1.tick_params(axis='y',which='both', right=True, labelright=True)
ax1.tick_params(axis='both', which='major', labelsize=10) 
ax1.set_title('Conserved quantities')
particle_mu, = ax1.plot([], [], color='red',label='Magnetic moment [J/T]')
particle_energy, = ax1.plot([], [], color='blue',label='Total energy [J]')
ax1.legend(loc=1)



# Downsample the trajectory to reduce the number of frames
mu1=mu[::step]
energy1=energy[::step]

# Initialization function for the animation
def init():
    particle_mu.set_data([], [])
    particle_energy.set_data([], [])
    return particle_mu, particle_energy

# Update function for the animation
def update(frame):
    # Update 2D plot
    particle_mu.set_data(time1[:frame], mu1[:frame])
    
    # Update 3D plot
    particle_energy.set_data(time1[:frame], energy1[:frame])
    
    return particle_mu, particle_energy

# Create the animation
ani = animation.FuncAnimation(
    fig, update, frames=len(x1), init_func=init, blit=False, interval=20)

# Save the combined animation to a file
writer = FFMpegWriter(fps=50, metadata=dict(artist='Me'), bitrate=1800)
ani.save('/Users/mv58/Documents/PhD/Plasma/Mirror25/conservation.mp4', writer=writer)
plt.show()

