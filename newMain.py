import matplotlib.pyplot as plt
# For PyCharm users; use TkAgg
#import matplotlib; matplotlib.use("TkAgg")
import numpy as np
from matplotlib import animation
from matplotlib.animation import FuncAnimation

# Orbital parameters of the satellite
AU = 149597870700 #meters in 1 AU
my_0 = 1.327E20 #heliocentric gravitational constant m^3/s^-2
my = my_0 / AU
a = 0.5959  # semimajor axis in AU
e = 0.5301  # eccentricity
inc = 33  # inclination in degrees
peri_arg = 0  # argument of periapsis (see link) in degrees
# https://mars.nasa.gov/mgs/status/nav/orbparam.html#:~:text=The%20argument%20of%20periapsis%20is,periapsis%20was%20about%20148.44%20degrees.
period = 168  # orbital period in days


frames = period * 2 # number of frames in the animation
time_res = 1 # stepsize

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
    # theta = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    beta = e / 1 + np.sqrt(1-e**2)
    theta = E + 2 * np.arctan2(beta * np.sin(E), 1 - beta * np.cos(E))

    # Distance between the focus of attraction and the orbiting body

    r = a * (1 - e ** 2) / (1 + e * np.cos(theta))


    # sin_v = np.sqrt(1 - e ** 2) * np.sin(E) / (1 - e * np.cos(E))
    # cos_v = (np.cos(E) - e) / (1 - e * np.cos(E))
    # # true anomaly?
    # v = np.arctan2(sin_v, cos_v)

    # Calculate the position vector in the orbital frame
    x_orb = r * np.cos(theta + np.radians(peri_arg))
    y_orb = r * np.sin(theta + np.radians(peri_arg))
    z_orb = 0

    # Convert to heliocentric coordinates
    x = x_orb
    y = y_orb * np.cos(np.radians(inc))
    z = y_orb * np.sin(np.radians(inc))

    # Velocity try 3
    v = np.sqrt(my * ((2/r) - (1/a))) # in m / s
    v_mag = v / 1000 # in km / s

    # Calculate the velocity of the satellite in Cartesian coordinates
    n = 2 * np.pi / (period * 24 * 60 * 60)
    vx = -r * n * np.sin(theta)
    vy = r * n * (e + np.cos(theta)) * np.cos(inc) / np.sqrt(1 - e ** 2)
    vz = r * n * (e + np.cos(theta)) * np.sin(inc) / np.sqrt(1 - e ** 2)
    # print('x:'+ str(vx) +'y:' + str(vy) + 'z:'+ str(vz))
    # v_mag = r * np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)

    return x, y, z, v_mag, r


# Create a 3D plot and add the Sun as a fixed point at the origin
# fig = plt.figure(figsize=(12, 6))
# ax = fig.add_subplot(121, projection='3d')
# ax2 = fig.add_subplot(222)
# ax3 = fig.add_subplot(224)
# ax.set_xlim3d(-1, 1)
# ax.set_ylim3d(-1, 1)
# ax.set_zlim3d(-1, 1)
# ax.scatter(0, 0, 0, color='yellow', s=200)
#
# legend = plt.legend(loc=2, prop={'size': 6})
#
# def init():
#     line.set_data([], [])
#     return line,


# Create the velocity subplot
#
# line, = ax2.plot([], [], 'r')
# xdata, ydata = [], []
#
#
# line2, = ax3.plot([], [], 'b')
# xdata2, ydata2 = [], []


# Define the update function that will be called for each frame
positions_x_y_z = [], [], []
velocities_x_y_z = [], [], []

# def update(frame):
#     positions = []
#     velocities = []
#     # Calculate the position of the satellite at the current time
#     t = frame * time_res
#     x, y, z, v_mag, r = calc_satellite_pos(t)
#
#     xdata.append(t)
#     ydata.append(v_mag)
#     ydata2.append(r)
#     # Should record the array of positions to animate the tail of the satellite
#     # TODO: does not work
#     positions.append((x, y, z))
#     velocities.append(v_mag)
#
#     positions_x_y_z[0].append(x)
#     positions_x_y_z[1].append(y)
#     positions_x_y_z[2].append(z)
#
#     velocities_x_y_z[0].append(t)
#     velocities_x_y_z[1].append(v_mag)
#     velocities_x_y_z[2].append(r)
#
#
#     vel = 'velocity = ' + str(round(v_mag, 2))
#     # if len(positions) > 1:
#     #     xs, ys, zs = zip(*positions)
#     #     ax.plot(xs, ys, zs, color='gray')
#
#     # Clear the previous frame and update the plot with the new position of the satellite
#     ax.clear()
#     ax.set_xlim3d(-1, 1)
#     ax.set_ylim3d(-1, 1)
#     ax.set_zlim3d(-1, 1)
#     ax.scatter(0, 0, 0, color='yellow', s=200, label='Sun')
#     ax.scatter(x, y, z, color='blue', s=10, label='Satellite')
#
#     ax.legend()
#     # ax.view_init(30, t * 360 / period)
#
#     # Update the velocity subplot
#     ax2.set_xlim(t-t, 5 + t)
#     ax2.set_ylim(5, 80)
#     print(v_mag)
#     ax2.set_xlabel('Time')
#     ax2.set_ylabel('Velocity')
#     line.set_data(xdata, ydata)
#
#     # Update the radius subplot
#     ax3.set_xlim(t-t, 5 + t)
#     ax3.set_ylim(0, 1)
#     print(v_mag)
#     ax3.set_xlabel('Time')
#     ax3.set_ylabel('Radius')
#     line2.set_data(xdata, ydata2)
#
#
#     # Define the XY plane
#     X, Y = np.meshgrid(np.arange(-1, 1, 0.1), np.arange(-1, 1, 0.1))
#     Z = np.zeros(X.shape)
#     ax.plot_surface(X, Y, Z, alpha=0.5)

# Create the animation using the FuncAnimation function
#anim = FuncAnimation(fig, update, init_func=init, frames=frames, interval=20)
# Display the animation
#plt.show()

# Exchange with the location of your desire
#f = r"c://Users/Georg/Documents/Uni/Space Sciences MA/Measurement Methods in Space Physics/animate_func.gif"
#writergif = animation.PillowWriter(fps=300 / 6)
#anim.save(f, writer=writergif)



def plot_2D_orbits():
    # make time array of integers for calculating the orbit
    x_pos_array = []
    y_pos_array = []
    z_pos_array = []

    # make integer array from 0 to 100
    time = np.arange(0, 170, 1)

    for i in range(len(time)):
        x, y, z, v_mag, r = calc_satellite_pos(i)
        x_pos_array.append(x)
        y_pos_array.append(y)
        z_pos_array.append(z)





    # plot the orbit in 3 subplots
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(10, 12))

    # title for all
    fig.suptitle('Orbit of the satellite split into x, y and z components for 1 orbital period', fontsize=16)

    ax1.plot(time, x_pos_array, color='blue')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('x position')

    ax2.plot(time, y_pos_array, color='green')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('y position')

    ax3.plot(time, z_pos_array, color='red')
    ax3.set_xlabel('Time')
    ax3.set_ylabel('z position')

    ax4.plot(time, x_pos_array, color='blue', label='x position')
    ax4.plot(time, y_pos_array, color='green', label='y position')
    ax4.plot(time, z_pos_array, color='red', label='z position')
    ax4.set_xlabel('Time')
    ax4.set_ylabel('Position in x, y, z')

    return


def draw_views():
    x_pos_array = []
    y_pos_array = []
    z_pos_array = []

    # make integer array from 0 to 100
    time = np.arange(0, 168*30, 1)

    for i in range(len(time)):
        x, y, z, v_mag, r = calc_satellite_pos(i)
        x_pos_array.append(x)
        y_pos_array.append(y)
        z_pos_array.append(z)

    # find first 168 values
    x_pos_array_1 = x_pos_array[0:168]
    y_pos_array_1 = y_pos_array[0:168]
    z_pos_array_1 = z_pos_array[0:168]

    #find last 168 values
    x_pos_array_2 = x_pos_array[-168:]
    y_pos_array_2 = y_pos_array[-168:]
    z_pos_array_2 = z_pos_array[-168:]

    # plot the orbit in 3 subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))
    ax1.plot(x_pos_array_1, y_pos_array_1, color='blue', label='First orbit')
    ax1.plot(x_pos_array_2, y_pos_array_2, color='red', linestyle='dashed', label='Last orbit')
    ax1.set_xlabel('x position')
    ax1.set_ylabel('y position')
    ax1.scatter(0, 0, color='yellow', s=200, label='Sun')
    ax1.legend()
    ax1.set_title('Top view of the orbit of the satellite')

    ax2.plot(x_pos_array_1, z_pos_array_1, color='blue', label='First orbit')
    ax2.plot(x_pos_array_2, z_pos_array_2, color='red', linestyle='dashed', label='Last orbit')
    ax2.set_xlabel('x position')
    ax2.set_ylabel('z position')
    ax2.scatter(0, 0, color='yellow', s=200, label='Sun')
    ax2.legend()
    ax2.set_title('Side view of the orbit of the satellite')

    ax3.plot(y_pos_array_1, z_pos_array_1, color='blue', label='First orbit')
    ax3.plot(y_pos_array_2, z_pos_array_2, color='red', linestyle='dashed', label='Last orbit')
    ax3.set_xlabel('y position')
    ax3.set_ylabel('z position')
    ax3.scatter(0, 0, color='yellow', s=200, label='Sun')
    ax3.legend()
    ax3.set_title('Front view of the orbit of the satellite')

    plt.show()
plot_2D_orbits()
draw_views()







