import matplotlib.pyplot as plt
# For PyCharm users; use TkAgg
import matplotlib; matplotlib.use("TkAgg")
import numpy as np
from matplotlib import animation
from matplotlib.animation import FuncAnimation

# Orbital parameters of the satellite
a = 0.5959  # semimajor axis in AU
e = 0.5301  # eccentricity
inc = 33  # inclination in degrees
peri_arg = 0  # argument of periapsis (see link) in degrees
# https://mars.nasa.gov/mgs/status/nav/orbparam.html#:~:text=The%20argument%20of%20periapsis%20is,periapsis%20was%20about%20148.44%20degrees.
period = 168  # orbital period in days


def calc_satellite_pos(t):
    # Convert time to radians
    M = 2 * np.pi * t / period

    # Solve Keplers equation for the eccentric anomaly
    E = M
    while True:
        E_new = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        if abs(E_new - E) < 1e-6:
            E = E_new
            break
        E = E_new

    # Calculate the true anomaly and distance
    theta = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    r = a * (1 - e ** 2) / (1 + e * np.cos(theta))
    sin_v = np.sqrt(1 - e ** 2) * np.sin(E) / (1 - e * np.cos(E))
    cos_v = (np.cos(E) - e) / (1 - e * np.cos(E))
    v = np.arctan2(sin_v, cos_v)

    # Calculate the position vector in the orbital frame
    x_orb = r * np.cos(theta + np.radians(peri_arg))
    y_orb = r * np.sin(theta + np.radians(peri_arg))
    z_orb = 0

    # Convert to heliocentric coordinates
    x = x_orb
    y = y_orb * np.cos(np.radians(inc))
    z = y_orb * np.sin(np.radians(inc))

    # Calculate the velocity of the satellite in Cartesian coordinates
    n = 2 * np.pi / (period * 24 * 60 * 60)
    vx = -r * n * np.sin(v)
    vy = r * n * (e + np.cos(v)) * np.cos(inc) / np.sqrt(1 - e ** 2)
    vz = r * n * (e + np.cos(v)) * np.sin(inc) / np.sqrt(1 - e ** 2)
    # print('x:'+ str(vx) +'y:' + str(vy) + 'z:'+ str(vz))
    v_mag = r * np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)

    return x, y, z, v_mag, r


# Create a 3D plot and add the Sun as a fixed point at the origin
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.set_xlim3d(-1, 1)
ax.set_ylim3d(-1, 1)
ax.set_zlim3d(-1, 1)
ax.scatter(0, 0, 0, color='yellow', s=200)

legend = plt.legend(loc=2, prop={'size': 6})

def init():
    line.set_data([], [])
    return line,


# Create the velocity subplot
ax2 = fig.add_subplot(1, 2, 2)
line, = ax2.plot([], [], 'r')
xdata, ydata = [], []

ax3 = fig.add_subplot(2, 2, 2)
line2, = ax3.plot([], [], 'b')
xdata2, ydata2 = [], []


# Define the update function that will be called for each frame
def update(frame):
    positions = []
    velocities = []
    # Calculate the position of the satellite at the current time
    t = frame * 5
    x, y, z, v_mag, r = calc_satellite_pos(t)

    xdata.append(t)
    ydata.append(v_mag)
    ydata2.append(r)
    # Should record the array of positions to animate the tail of the satellite
    # TODO: does not work
    positions.append((x, y, z))
    velocities.append(v_mag)

    vel = 'velocity = ' + str(round(v_mag, 2))
    # if len(positions) > 1:
    #     xs, ys, zs = zip(*positions)
    #     ax.plot(xs, ys, zs, color='gray')

    # Clear the previous frame and update the plot with the new position of the satellite
    ax.clear()
    ax.set_xlim3d(-1, 1)
    ax.set_ylim3d(-1, 1)
    ax.set_zlim3d(-1, 1)
    ax.scatter(0, 0, 0, color='yellow', s=200, label='Sun')
    ax.scatter(x, y, z, color='blue', s=10, label='Satellite')
    ax.axvline(x=0,ymax=1, ymin=-1, color='black')

    ax.legend()
    # ax.view_init(30, t * 360 / period)

    # Update the velocity subplot
    ax2.set_xlim(t-t, 5 + t)
    ax2.set_ylim(0, 1E-6)
    print(v_mag)
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Velocity')
    line.set_data(xdata, ydata)

    # Update the radius subplot
    ax3.set_xlim(t-t, 5 + t)
    ax3.set_ylim(0, 1)
    print(v_mag)
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Radius')
    line2.set_data(xdata, ydata2)


    # Define the XY plane
    X, Y = np.meshgrid(np.arange(-1, 1, 0.1), np.arange(-1, 1, 0.1))
    Z = np.zeros(X.shape)
    ax.plot_surface(X, Y, Z, alpha=0.5)

# Create the animation using the FuncAnimation function
anim = FuncAnimation(fig, update, init_func=init, frames=100, interval=20)
# Display the animation
plt.show()
# Exchange with the location of your desire
f = r"c://Users/Georg/Documents/Uni/Space Sciences MA/Measurement Methods in Space Physics/animate_func.gif"
writergif = animation.PillowWriter(fps=30 / 6)
anim.save(f, writer=writergif)
