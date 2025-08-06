#!/Users/subhodeep/venv/bin/python
# import pandas as pd
# import matplotlib.pyplot as plt
# import glob

# files = glob.glob("*.csv")

# fig, axs = plt.subplots(len(files), 1, figsize=(8, 3*len(files)))

# for ax, file in zip(axs, files):
#     data = pd.read_csv(file)
#     ax.plot(data.x, data.y, '-o', markersize=2)
#     ax.set_title(file)
#     ax.set_xlabel("x (or s)")
#     ax.set_ylabel("y (or d)")
#     ax.grid(True)

# plt.tight_layout()
# plt.show()




# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation
# from matplotlib.patches import Circle

# def load(fname): return pd.read_csv(fname).values

# # Load Existing Data
# path = load("mule_nc.csv")
# mule_nc, joint_nc, trailer_nc = load("mule_nc.csv"), load("joint_nc.csv"), load("trailer_nc.csv")
# mule_mpc, joint_mpc, trailer_mpc = load("mule_mpc.csv"), load("joint_mpc.csv"), load("trailer_mpc.csv")
# mule_rev, joint_rev, trailer_rev = load("mule_reverse.csv"), load("joint_reverse.csv"), load("trailer_reverse.csv")
# mule_obs, joint_obs, trailer_obs = load("mule_obstacles.csv"), load("joint_obstacles.csv"), load("trailer_obstacles.csv")

# # ✅ Load New Safe Path and Trajectories
# safe_path = load("safe_path.csv")
# mule_avoid, joint_avoid, trailer_avoid = load("mule_avoid.csv"), load("joint_avoid.csv"), load("trailer_avoid.csv")

# # Obstacles (unchanged)
# obstacles = [(10.0, 0.0, 0.5), (-10.0, 0.0, 0.8), (5.0, 4.0, 0.8)]

# def animate_case(path_points, mule, joint, trailer, title, draw_obs=False):
#     fig, ax = plt.subplots()
#     ax.plot(path_points[:,0], path_points[:,1], 'b--', label="Reference Path")

#     # Points
#     mule_point, = ax.plot([], [], 'ro', label="Mule")
#     joint_point, = ax.plot([], [], 'yo', label="Joint")
#     trailer_point, = ax.plot([], [], 'ko', label="Trailer")

#     # Lines for joint visualization
#     line1, = ax.plot([], [], 'r-', lw=1.5)  # Mule → Joint
#     line2, = ax.plot([], [], 'k-', lw=1.5)  # Joint → Trailer

#     ax.axis('equal'); ax.legend(); ax.set_title(title)

#     obs_patches = []
#     if draw_obs:
#         for ox, oy, r in obstacles:
#             c = Circle((ox, oy), r, color='r', alpha=0.3)
#             ax.add_patch(c)
#             obs_patches.append(c)

#     def init():
#         mule_point.set_data([], [])
#         joint_point.set_data([], [])
#         trailer_point.set_data([], [])
#         line1.set_data([], [])
#         line2.set_data([], [])
#         return mule_point, joint_point, trailer_point, line1, line2, *obs_patches

#     def update(i):
#         mule_point.set_data([mule[i,0]], [mule[i,1]])
#         joint_point.set_data([joint[i,0]], [joint[i,1]])
#         trailer_point.set_data([trailer[i,0]], [trailer[i,1]])

#         line1.set_data([mule[i,0], joint[i,0]], [mule[i,1], joint[i,1]])
#         line2.set_data([joint[i,0], trailer[i,0]], [joint[i,1], trailer[i,1]])

#         return mule_point, joint_point, trailer_point, line1, line2, *obs_patches

#     ani = animation.FuncAnimation(fig, update, frames=len(mule),
#                                   init_func=init, blit=True, interval=30)
#     plt.show()

# # ---- Existing Static Plot for Old Obstacle Case ----
# plt.plot(mule_obs[:,0], mule_obs[:,1], 'r-', label='Mule Obs')
# plt.plot(joint_obs[:,0], joint_obs[:,1], 'y-', label='Joint Obs')
# plt.plot(trailer_obs[:,0], trailer_obs[:,1], 'k-', label='Trailer Obs')
# plt.legend(); plt.axis('equal'); plt.show()

# # ---- Existing Animations ----
# animate_case(path, mule_nc, joint_nc, trailer_nc, "Forward (No Controller)")
# animate_case(path, mule_mpc, joint_mpc, trailer_mpc, "Forward (MPC Controller)")
# animate_case(path, mule_rev, joint_rev, trailer_rev, "Reverse (Mule Pushing Trailer)")
# animate_case(path, mule_obs, joint_obs, trailer_obs, "Forward with Old Obstacle Avoidance", draw_obs=True)

# # ✅ ---- New Animation for Precomputed Safe Path ----
# animate_case(safe_path, mule_avoid, joint_avoid, trailer_avoid, "Forward with Precomputed Safe Path", draw_obs=True)





import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle

def load(fname): return pd.read_csv(fname).values

# Load Data
# ✅ Load actual reference path
path = load("reference_path.csv")
mule_nc, joint_nc, trailer_nc = load("mule_nc.csv"), load("joint_nc.csv"), load("trailer_nc.csv")
mule_mpc, joint_mpc, trailer_mpc = load("mule_mpc.csv"), load("joint_mpc.csv"), load("trailer_mpc.csv")
mule_rev, joint_rev, trailer_rev = load("mule_reverse.csv"), load("joint_reverse.csv"), load("trailer_reverse.csv")

# ✅ Updated Obstacle Avoidance Data (use only *_avoid.csv)
safe_path = load("safe_path.csv")
mule_avoid, joint_avoid, trailer_avoid = load("mule_avoid.csv"), load("joint_avoid.csv"), load("trailer_avoid.csv")

# Obstacles
obstacles = [(10.0, 0.0, 0.5), (-10.0, 0.0, 0.8), (5.0, 4.0, 0.8)]

def animate_case(path_points, mule, joint, trailer, title, draw_obs=False):
    fig, ax = plt.subplots()
    ax.plot(path_points[:,0], path_points[:,1], 'b--', label="Reference Path")

    mule_point, = ax.plot([], [], 'ro', label="Mule")
    joint_point, = ax.plot([], [], 'yo', label="Joint")
    trailer_point, = ax.plot([], [], 'ko', label="Trailer")

    line1, = ax.plot([], [], 'r-', lw=1.5)  # Mule → Joint
    line2, = ax.plot([], [], 'k-', lw=1.5)  # Joint → Trailer

    ax.axis('equal'); ax.legend(); ax.set_title(title)

    obs_patches = []
    if draw_obs:
        for ox, oy, r in obstacles:
            c = Circle((ox, oy), r, color='r', alpha=0.3)
            ax.add_patch(c)
            obs_patches.append(c)

    def init():
        mule_point.set_data([], [])
        joint_point.set_data([], [])
        trailer_point.set_data([], [])
        line1.set_data([], [])
        line2.set_data([], [])
        return mule_point, joint_point, trailer_point, line1, line2, *obs_patches

    def update(i):
        mule_point.set_data([mule[i,0]], [mule[i,1]])
        joint_point.set_data([joint[i,0]], [joint[i,1]])
        trailer_point.set_data([trailer[i,0]], [trailer[i,1]])

        line1.set_data([mule[i,0], joint[i,0]], [mule[i,1], joint[i,1]])
        line2.set_data([joint[i,0], trailer[i,0]], [joint[i,1], trailer[i,1]])

        return mule_point, joint_point, trailer_point, line1, line2, *obs_patches

    ani = animation.FuncAnimation(fig, update, frames=len(mule),
                                  init_func=init, blit=True, interval=30)
    plt.show()

# ---- Static Plot for Avoidance Path ----
plt.plot(mule_avoid[:,0], mule_avoid[:,1], 'r-', label='Mule Avoid')
plt.plot(joint_avoid[:,0], joint_avoid[:,1], 'y-', label='Joint Avoid')
plt.plot(trailer_avoid[:,0], trailer_avoid[:,1], 'k-', label='Trailer Avoid')
plt.legend(); plt.axis('equal'); plt.show()

# ---- Animations ----
animate_case(path, mule_nc, joint_nc, trailer_nc, "Forward (No Controller)")
animate_case(path, mule_mpc, joint_mpc, trailer_mpc, "Forward (MPC Controller)")
animate_case(path, mule_rev, joint_rev, trailer_rev, "Reverse (Mule Pushing Trailer)")
animate_case(safe_path, mule_avoid, joint_avoid, trailer_avoid, "Forward with Precomputed Safe Path", draw_obs=True)
