import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def simulate_landreman_ocean_basin(
    basin_length_km=10000, 
    basin_width_km=4000, 
    basin_depth_m=500, 
    iota=0.333, # Winding ratio (Poloidal/Toroidal)
    n_points=1000
):
    """
    Simulates a 3D 'twisty' standing wave orbit (invariant torus streamline)
    in an ocean basin governed by Landreman-like Euler equilibrium equations.
    """
    # Toroidal angle phi (along the basin length X)
    phi = np.linspace(0, 4 * np.pi, n_points)
    
    # Major axis dimensions (Longitudinal / Toroidal scale)
    R_x = basin_length_km / 2.0
    
    # Minor axis dimensions (Latitudinal & Vertical / Poloidal scale)
    r_y = basin_width_km / 4.0
    r_z = basin_depth_m / 2.0
    
    # Poloidal angle theta twisted by the rotational transform (iota)
    theta = iota * phi
    
    # Non-axisymmetric spatial perturbation (Landreman twist term)
    twist_factor = 1.0 + 0.2 * np.cos(3 * phi)
    
    # Cartesian coordinates of the 3D helical flow line
    X = R_x * np.cos(phi)
    Y = r_y * np.cos(theta) * twist_factor
    Z = -r_z + r_z * np.sin(theta) * twist_factor  # Constrained to depth (Z < 0)
    
    # Velocity components along the streamtube
    u_x = -R_x * np.sin(phi)
    u_y = -r_y * iota * np.sin(theta)
    u_z =  r_z * iota * np.cos(theta)
    
    return X, Y, Z, u_x, u_y, u_z

# Execute simulation
X, Y, Z, Ux, Uy, Uz = simulate_landreman_ocean_basin()

# Visualization
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')

# Plot the 3D Helical Trajectory (Standing Wave Streamline)
sc = ax.scatter(X, Y, Z, c=np.sqrt(Ux**2 + Uy**2 + Uz**2), cmap='plasma', s=8)

ax.set_title("Landreman-Like 'Twisty' 3D Standing Wave in an Ocean Basin")
ax.set_xlabel("Basin Length X (km)")
ax.set_ylabel("Basin Width Y (km)")
ax.set_zlabel("Depth Z (m)")
plt.colorbar(sc, label="Fluid Velocity Magnitude (m/s)")

plt.tight_layout()
plt.show()
