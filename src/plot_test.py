#!/Users/subhodeep/venv/bin/python

# import pandas as pd
# import matplotlib.pyplot as plt

# # Load data
# path = pd.read_csv("mule_nc.csv")
# trailer_nc = pd.read_csv("trailer_nc.csv")
# trailer_mpc = pd.read_csv("trailer_mpc.csv")
# err_nc = pd.read_csv("err_nc.csv")
# err_mpc = pd.read_csv("err_mpc.csv")

# fig, ax = plt.subplots(1,2,figsize=(12,5))
# ax[0].plot(path.x, path.y, 'b--', label="Path")
# ax[0].plot(trailer_nc.x, trailer_nc.y, 'r', label="Trailer NC")
# ax[0].plot(trailer_mpc.x, trailer_mpc.y, 'g', label="Trailer MPC")
# ax[0].set_title("Trajectories"); ax[0].legend(); ax[0].axis('equal')

# ax[1].plot(err_nc.error, 'r', label="Error NC")
# ax[1].plot(err_mpc.error, 'g', label="Error MPC")
# ax[1].set_title("Trailer Error"); ax[1].legend()

# plt.show()



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# -------- Load Data --------
path = pd.read_csv("mule_nc.csv")  # using mule_nc just to get reference path shape
path_points = path[['x','y']].values

mule_nc = pd.read_csv("mule_nc.csv").values
joint_nc = pd.read_csv("joint_nc.csv").values
trailer_nc = pd.read_csv("trailer_nc.csv").values

mule_mpc = pd.read_csv("mule_mpc.csv").values
joint_mpc = pd.read_csv("joint_mpc.csv").values
trailer_mpc = pd.read_csv("trailer_mpc.csv").values

# -------- Animation Function --------
def animate_three_points(path_points, mule_traj, joint_traj, trailer_traj, title):
    fig, ax = plt.subplots()
    ax.plot(path_points[:,0], path_points[:,1], 'b--', label="Figure-8 Path")
    mule_point, = ax.plot([], [], 'ro', label="Mule")
    joint_point, = ax.plot([], [], 'yo', label="Joint")
    trailer_point, = ax.plot([], [], 'ko', label="Trailer")
    
    ax.axis('equal'); ax.legend(); ax.set_title(title)

    def init():
        mule_point.set_data([], [])
        joint_point.set_data([], [])
        trailer_point.set_data([], [])
        return mule_point, joint_point, trailer_point

    def update(i):
        mule_point.set_data([mule_traj[i,0]], [mule_traj[i,1]])
        joint_point.set_data([joint_traj[i,0]], [joint_traj[i,1]])
        trailer_point.set_data([trailer_traj[i,0]], [trailer_traj[i,1]])
        return mule_point, joint_point, trailer_point


    ani = animation.FuncAnimation(fig, update, frames=len(mule_traj), init_func=init, blit=True, interval=30)
    plt.show()

# -------- Run Animations --------
animate_three_points(path_points, mule_nc, joint_nc, trailer_nc, "Without MPC Controller")
animate_three_points(path_points, mule_mpc, joint_mpc, trailer_mpc, "With MPC Controller")
