#!/Users/subhodeep/venv/bin/python
# import numpy as np
# import matplotlib.pyplot as plt

# # ================== Frenet Frame Generator ==================
# def generate_frenet_frame(parametric_func, t_range=(0, 2*np.pi), num_points=500):
#     """
#     Generate Frenet coordinates (s, d) for a given parametric curve defined by parametric_func(t) -> (x, y).
#     Plot both the curve in Cartesian space with Frenet frames and its transformed version in Frenet space.
#     """
#     t = np.linspace(t_range[0], t_range[1], num_points)
#     points = np.array([parametric_func(tt) for tt in t])
    
#     # Compute tangent vectors
#     dt = t[1] - t[0]
#     dx = np.gradient(points[:,0], dt)
#     dy = np.gradient(points[:,1], dt)
#     tangents = np.vstack((dx, dy)).T
#     tangent_norms = np.linalg.norm(tangents, axis=1).reshape(-1,1)
#     tangents_unit = tangents / tangent_norms
    
#     # Compute normals by rotating tangents by -90 degrees
#     normals_unit = np.zeros_like(tangents_unit)
#     normals_unit[:,0] = -tangents_unit[:,1]
#     normals_unit[:,1] = tangents_unit[:,0]
    
#     # Compute s (arc-length) and d (lateral offset = 0 for reference curve)
#     s = np.zeros(num_points)
#     for i in range(1, num_points):
#         s[i] = s[i-1] + np.linalg.norm(points[i] - points[i-1])
#     d = np.zeros(num_points)
#     frenet_coords = np.vstack((s, d)).T
    
#     # ---------- Plot 1: Cartesian space with Frenet frames ----------
#     plt.figure(figsize=(8,6))
#     plt.plot(points[:,0], points[:,1], 'b', label='Curve in Cartesian')
#     for i in range(0, num_points, num_points//20):
#         p = points[i]
#         t_vec = tangents_unit[i] * 0.5
#         n_vec = normals_unit[i] * 0.5
#         plt.arrow(p[0], p[1], t_vec[0], t_vec[1], color='g', head_width=0.1)
#         plt.arrow(p[0], p[1], n_vec[0], n_vec[1], color='r', head_width=0.1)
#     plt.axis('equal')
#     plt.legend(["Curve", "Tangent (green)", "Normal (red)"])
#     plt.title("Frenet Frame Visualization in Cartesian Space")
#     # plt.show()
    
#     # ---------- Plot 2: Frenet space (s vs d) ----------
#     plt.figure(figsize=(8,4))
#     plt.plot(frenet_coords[:,0], frenet_coords[:,1], 'm')
#     plt.title("Curve Representation in Frenet Frame (s-d space)")
#     plt.xlabel("s (arc-length)")
#     plt.ylabel("d (lateral offset)")
#     plt.grid(True)
#     plt.show()
    
#     return frenet_coords

# # ================== Example Usage ==================
# # Parametric curve: figure-8 using sin and cos
# def figure_eight(t):
#     return np.array([10*np.sin(t), 10*np.sin(t)*np.cos(t)])

# frenet_coords = generate_frenet_frame(figure_eight)
# print("Frenet coordinates computed:", frenet_coords.shape)



# ================================================================================================================================
# ================================================================================================================================
# ================================================================================================================================
# ================================================================================================================================



# import numpy as np
# import matplotlib.pyplot as plt

# # ================== Frenet Frame Generator ==================
# def generate_frenet_frame(parametric_func, t_range=(0, 2*np.pi), num_points=500):
#     t = np.linspace(t_range[0], t_range[1], num_points)
#     points = np.array([parametric_func(tt) for tt in t])
    
#     dt = t[1] - t[0]
#     dx = np.gradient(points[:,0], dt)
#     dy = np.gradient(points[:,1], dt)
#     tangents = np.vstack((dx, dy)).T
#     tangent_norms = np.linalg.norm(tangents, axis=1).reshape(-1,1)
#     tangents_unit = tangents / tangent_norms

#     normals_unit = np.zeros_like(tangents_unit)
#     normals_unit[:,0] = -tangents_unit[:,1]
#     normals_unit[:,1] = tangents_unit[:,0]

#     s = np.zeros(num_points)
#     for i in range(1, num_points):
#         s[i] = s[i-1] + np.linalg.norm(points[i] - points[i-1])
#     d = np.zeros(num_points)
#     frenet_coords = np.vstack((s, d)).T

#     # Plot Cartesian
#     plt.figure(figsize=(8,6))
#     plt.plot(points[:,0], points[:,1], 'b', label='Curve in Cartesian')
#     plt.axis('equal'); plt.title("Frenet Frame in Cartesian"); plt.show()

#     # Plot Frenet
#     plt.figure(figsize=(8,4))
#     plt.plot(frenet_coords[:,0], frenet_coords[:,1], 'm')
#     plt.title("Curve in Frenet Space (s-d)"); plt.xlabel("s"); plt.ylabel("d"); plt.grid(); plt.show()
#     return points, frenet_coords

# # ================== Simple Pure Pursuit Controller ==================
# def pure_pursuit(current_pos, path, lookahead=1.0, v=0.2):
#     dists = np.linalg.norm(path - current_pos, axis=1)
#     idx = np.argmin(dists)
#     goal_idx = min(idx + 5, len(path)-1)
#     goal = path[goal_idx]
#     direction = (goal - current_pos)/np.linalg.norm(goal - current_pos)
#     return current_pos + v * direction

# # ================== Obstacle Avoidance ==================
# def add_obstacle(path, obstacle_pos, avoidance_offset=3.0):
#     new_path = []
#     for p in path:
#         if np.linalg.norm(p - obstacle_pos) < 5:
#             new_path.append(p + [0, avoidance_offset])
#         else:
#             new_path.append(p)
#     return np.array(new_path)

# def simulate_tracking(path):
#     pos = path[0].copy(); traj = [pos]
#     for _ in range(1000):
#         pos = pure_pursuit(pos, path)
#         traj.append(pos)
#         if np.linalg.norm(pos - path[-1]) < 0.5: break
#     return np.array(traj)

# # ================== Trailer Kinematics ==================
# def trailer_kinematics(pos, theta, phi, v=0.2, L=2.0, D=2.0, dt=0.1):
#     x, y = pos
#     x += v*np.cos(theta)*dt
#     y += v*np.sin(theta)*dt
#     theta += (v/L)*np.tan(phi)*dt
#     return np.array([x, y]), theta

# def simulate_trailer(path):
#     pos = path[0].copy(); theta = 0.0; phi = 0.0; traj = [pos]
#     for _ in range(1000):
#         pos, theta = trailer_kinematics(pos, theta, phi)
#         traj.append(pos)
#         if np.linalg.norm(pos - path[-1]) < 0.5: break
#     return np.array(traj)

# # ================== Example Usage ==================
# def figure_eight(t):
#     return np.array([10*np.sin(t), 10*np.sin(t)*np.cos(t)])

# # Generate original path
# path_cart, frenet_coords = generate_frenet_frame(figure_eight)

# # --- Part 3a: Obstacle Avoidance ---
# obstacle = np.array([0,0])
# path_avoid = add_obstacle(path_cart, obstacle)
# traj_avoid = simulate_tracking(path_avoid)
# plt.figure(); plt.plot(path_cart[:,0], path_cart[:,1],'b',label='Original');
# plt.plot(path_avoid[:,0], path_avoid[:,1],'c--',label='Avoidance Path');
# plt.plot(traj_avoid[:,0], traj_avoid[:,1],'r',label='Robot Trajectory');
# plt.scatter(*obstacle,color='k',label='Obstacle');
# plt.legend(); plt.axis('equal'); plt.title("Obstacle Avoidance"); plt.show()

# # --- Part 3b: Trailer Tracking in Frenet Frame ---
# traj_trailer = simulate_trailer(path_cart)
# plt.figure(); plt.plot(path_cart[:,0], path_cart[:,1],'b',label='Reference');
# plt.plot(traj_trailer[:,0], traj_trailer[:,1],'g',label='Trailer Trajectory');
# plt.axis('equal'); plt.legend(); plt.title("Trailer Tracking Figure-8"); plt.show()

# # Plot Trailer Path in Frenet Space
# s_trailer = np.cumsum(np.r_[0, np.linalg.norm(np.diff(traj_trailer, axis=0), axis=1)])
# d_trailer = np.zeros_like(s_trailer)
# plt.figure(); plt.plot(s_trailer, d_trailer,'g'); plt.title("Trailer in Frenet Space (s-d)"); plt.xlabel("s"); plt.ylabel("d"); plt.grid(); plt.show()

# print("Obstacle avoidance and trailer tracking simulations complete.")



# ================================================================================================================================
# ================================================================================================================================
# ================================================================================================================================
# ================================================================================================================================



# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation

# # ================== Frenet Frame Generator ==================
# def generate_frenet_frame(parametric_func, t_range=(0, 2*np.pi), num_points=500):
#     t = np.linspace(t_range[0], t_range[1], num_points)
#     points = np.array([parametric_func(tt) for tt in t])
#     return points

# # ================== Simulate Mule, Joint, Trailer ==================
# def simulate_three_point_body(path, L=1.5):
#     mule_positions = []
#     joint_positions = []
#     trailer_positions = []
#     for i in range(len(path)):
#         mule = path[i]
#         if i == 0:
#             prev = path[i]
#         else:
#             prev = path[i-1]
#         direction = mule - prev
#         if np.linalg.norm(direction) < 1e-6:
#             direction = np.array([1.0,0.0])
#         else:
#             direction = direction / np.linalg.norm(direction)
#         # joint is L behind mule, trailer is L behind joint
#         joint = mule - L*direction
#         trailer = joint - L*direction
#         mule_positions.append(mule)
#         joint_positions.append(joint)
#         trailer_positions.append(trailer)
#     return np.array(mule_positions), np.array(joint_positions), np.array(trailer_positions)

# # ================== Animation ==================
# def animate_three_points(path, mule_traj, joint_traj, trailer_traj):
#     fig, ax = plt.subplots()
#     ax.plot(path[:,0], path[:,1], 'b--', label='Figure-8 Path')
#     mule_point, = ax.plot([], [], 'ro', label='Mule')
#     joint_point, = ax.plot([], [], 'yo', label='Joint')
#     trailer_point, = ax.plot([], [], 'ko', label='Trailer')
#     ax.legend(); ax.axis('equal'); ax.set_title("3-Point Body Tracking Figure-8")

#     def init():
#         mule_point.set_data([], [])
#         joint_point.set_data([], [])
#         trailer_point.set_data([], [])
#         return mule_point, joint_point, trailer_point

#     def update(i):
#         mule_point.set_data([mule_traj[i,0]], [mule_traj[i,1]])
#         joint_point.set_data([joint_traj[i,0]], [joint_traj[i,1]])
#         trailer_point.set_data([trailer_traj[i,0]], [trailer_traj[i,1]])
#         return mule_point, joint_point, trailer_point

#     ani = animation.FuncAnimation(fig, update, frames=len(mule_traj), init_func=init,
#                                   interval=20, blit=True, repeat=False)
#     plt.show()

# # ================== Example Usage ==================
# def figure_eight(t):
#     return np.array([10*np.sin(t), 10*np.sin(t)*np.cos(t)])

# path_cart = generate_frenet_frame(figure_eight)
# mule_traj, joint_traj, trailer_traj = simulate_three_point_body(path_cart, L=1.5)
# animate_three_points(path_cart, mule_traj, joint_traj, trailer_traj)




# ================================================================================================================================
# ================================================================================================================================
# ================================================================================================================================
# ================================================================================================================================




# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation

# # ================== Frenet Frame Generator ==================
# def generate_frenet_frame(parametric_func, t_range=(0, 2*np.pi), num_points=500):
#     t = np.linspace(t_range[0], t_range[1], num_points)
#     points = np.array([parametric_func(tt) for tt in t])
#     return points

# # ================== Simulate Mule, Flexible Joint, Trailer ==================
# def simulate_three_point_body(path, L=1.5, joint_flex=0.3):
#     mule_positions = []
#     joint_positions = []
#     trailer_positions = []
#     trailer_dir = np.array([-1.0, 0.0])  # initial trailer direction

#     for i in range(len(path)):
#         mule = path[i]
#         prev = path[i-1] if i > 0 else mule
#         mule_dir = mule - prev
#         if np.linalg.norm(mule_dir) < 1e-6:
#             mule_dir = np.array([1.0,0.0])
#         else:
#             mule_dir /= np.linalg.norm(mule_dir)

#         # Joint is not strictly behind mule; it can deviate slightly by an angle
#         angle_offset = joint_flex * np.sin(i * 0.1)  # oscillates for flexibility
#         rot = np.array([[np.cos(angle_offset), -np.sin(angle_offset)],
#                         [np.sin(angle_offset),  np.cos(angle_offset)]])
#         joint_dir = rot @ (-mule_dir)
#         joint = mule + L * joint_dir

#         # Trailer is behind joint along its own direction, which slowly aligns with joint_dir
#         trailer_dir = (0.9 * trailer_dir + 0.1 * joint_dir)
#         trailer_dir /= np.linalg.norm(trailer_dir)
#         trailer = joint + L * trailer_dir

#         mule_positions.append(mule)
#         joint_positions.append(joint)
#         trailer_positions.append(trailer)

#     return np.array(mule_positions), np.array(joint_positions), np.array(trailer_positions)

# # ================== Animation ==================
# def animate_three_points(path, mule_traj, joint_traj, trailer_traj):
#     fig, ax = plt.subplots()
#     ax.plot(path[:,0], path[:,1], 'b--', label='Figure-8 Path')
#     mule_point, = ax.plot([], [], 'ro', label='Mule')
#     joint_point, = ax.plot([], [], 'yo', label='Joint')
#     trailer_point, = ax.plot([], [], 'ko', label='Trailer')
#     ax.legend(); ax.axis('equal'); ax.set_title("Flexible Mule-Joint-Trailer Tracking Figure-8")

#     def init():
#         mule_point.set_data([], [])
#         joint_point.set_data([], [])
#         trailer_point.set_data([], [])
#         return mule_point, joint_point, trailer_point

#     def update(i):
#         mule_point.set_data([mule_traj[i,0]], [mule_traj[i,1]])
#         joint_point.set_data([joint_traj[i,0]], [joint_traj[i,1]])
#         trailer_point.set_data([trailer_traj[i,0]], [trailer_traj[i,1]])
#         return mule_point, joint_point, trailer_point

#     ani = animation.FuncAnimation(fig, update, frames=len(mule_traj), init_func=init,
#                                   interval=20, blit=True, repeat=False)
#     plt.show()

# # ================== Example Usage ==================
# def figure_eight(t):
#     return np.array([10*np.sin(t), 10*np.sin(t)*np.cos(t)])

# path_cart = generate_frenet_frame(figure_eight)
# mule_traj, joint_traj, trailer_traj = simulate_three_point_body(path_cart, L=1.5, joint_flex=0.4)
# animate_three_points(path_cart, mule_traj, joint_traj, trailer_traj)




# ================================================================================================================================
# ================================================================================================================================
# ================================================================================================================================
# ================================================================================================================================




# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation

# def generate_frenet_frame(parametric_func, t_range=(0, 2*np.pi), num_points=600):
#     t = np.linspace(t_range[0], t_range[1], num_points)
#     return np.array([parametric_func(tt) for tt in t])

# def mpc_controller(path, L=1.5, dt=0.05, gain=0.3):
#     mule_traj, joint_traj, trailer_traj = [], [], []
#     mule, trailer = path[0].copy(), path[0].copy() - np.array([2*L, 0])
#     joint = mule - L * np.array([1.0, 0.0])
#     trailer_dir = np.array([-1.0, 0.0])
#     for i in range(len(path)):
#         target = path[min(i+5, len(path)-1)]
#         mule_dir = target - mule
#         mule_dir /= np.linalg.norm(mule_dir)
#         mule += mule_dir * gain
#         dir_joint = mule - joint
#         dir_joint /= np.linalg.norm(dir_joint)
#         joint = mule - L * dir_joint
#         trailer_dir = 0.8 * trailer_dir + 0.2 * dir_joint
#         trailer_dir /= np.linalg.norm(trailer_dir)
#         trailer = joint - L * trailer_dir
#         mule_traj.append(mule.copy()); joint_traj.append(joint.copy()); trailer_traj.append(trailer.copy())
#     return np.array(mule_traj), np.array(joint_traj), np.array(trailer_traj)

# def simulate_three_point_body(path, L=1.5, joint_flex=0.3):
#     mule_traj, joint_traj, trailer_traj = [], [], []
#     trailer_dir = np.array([-1.0, 0.0])
#     for i, mule in enumerate(path):
#         prev = path[i-1] if i > 0 else mule
#         mule_dir = mule - prev
#         mule_dir = mule_dir / np.linalg.norm(mule_dir) if np.linalg.norm(mule_dir) > 1e-6 else np.array([1.0, 0.0])
#         angle_offset = joint_flex * np.sin(i * 0.05)
#         rot = np.array([[np.cos(angle_offset), -np.sin(angle_offset)], [np.sin(angle_offset), np.cos(angle_offset)]])
#         joint_dir = rot @ (-mule_dir)
#         joint = mule + L * joint_dir
#         trailer_dir = 0.9 * trailer_dir + 0.1 * joint_dir
#         trailer_dir /= np.linalg.norm(trailer_dir)
#         trailer = joint + L * trailer_dir
#         mule_traj.append(mule); joint_traj.append(joint); trailer_traj.append(trailer)
#     return np.array(mule_traj), np.array(joint_traj), np.array(trailer_traj)

# def compute_error(trailer_traj, path):
#     return np.array([np.min(np.linalg.norm(path - p, axis=1)) for p in trailer_traj])

# def animate_three_points(path, mule_traj, joint_traj, trailer_traj, title):
#     fig, ax = plt.subplots()
#     ax.plot(path[:, 0], path[:, 1], 'b--', label='Figure-8 Path')
#     mule_point, = ax.plot([], [], 'ro', label='Mule')
#     joint_point, = ax.plot([], [], 'yo', label='Joint')
#     trailer_point, = ax.plot([], [], 'ko', label='Trailer')
#     ax.legend(); ax.axis('equal'); ax.set_title(title)

#     def init():
#         mule_point.set_data([], []); joint_point.set_data([], []); trailer_point.set_data([], [])
#         return mule_point, joint_point, trailer_point

#     def update(i):
#         mule_point.set_data([mule_traj[i, 0]], [mule_traj[i, 1]])
#         joint_point.set_data([joint_traj[i, 0]], [joint_traj[i, 1]])
#         trailer_point.set_data([trailer_traj[i, 0]], [trailer_traj[i, 1]])
#         return mule_point, joint_point, trailer_point

#     ani = animation.FuncAnimation(fig, update, frames=len(mule_traj), init_func=init, interval=30, blit=True)
#     plt.show()

# def figure_eight(t):
#     return np.array([10 * np.sin(t), 10 * np.sin(t) * np.cos(t)])

# path = generate_frenet_frame(figure_eight)
# mule_nc, joint_nc, trailer_nc = simulate_three_point_body(path)
# mule_mpc, joint_mpc, trailer_mpc = mpc_controller(path)
# animate_three_points(path, mule_nc, joint_nc, trailer_nc, "Without MPC Controller")
# animate_three_points(path, mule_mpc, joint_mpc, trailer_mpc, "With MPC Controller")
# err_nc, err_mpc = compute_error(trailer_nc, path), compute_error(trailer_mpc, path)
# plt.figure(); plt.plot(err_nc, 'r', label='No Controller'); plt.plot(err_mpc, 'g', label='MPC Controller')
# plt.xlabel('Step'); plt.ylabel('Trailer Error'); plt.legend(); plt.title('Trailer Error Comparison'); plt.show()




# ================================================================================================================================
# ================================================================================================================================
# ================================================================================================================================
# ================================================================================================================================



import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def generate_frenet_frame(parametric_func, t_range=(0, 2*np.pi), num_points=600):
    t = np.linspace(t_range[0], t_range[1], num_points)
    return np.array([parametric_func(tt) for tt in t])

def mpc_controller(path, L=1.5, dt=0.05, gain=0.25, smooth=0.8):
    mule_traj, joint_traj, trailer_traj = [], [], []
    mule, trailer = path[0].copy(), path[0].copy() - np.array([2*L, 0])
    joint = mule - L * np.array([1.0, 0.0])
    trailer_dir = np.array([-1.0, 0.0])
    mule_vel = np.zeros(2)
    for i in range(len(path)):
        target = path[min(i+5, len(path)-1)]
        desired_dir = target - mule
        desired_dir /= np.linalg.norm(desired_dir)
        # Low-pass filter on velocity to smooth motion
        mule_vel = smooth * mule_vel + (1 - smooth) * (desired_dir * gain)
        mule += mule_vel
        dir_joint = mule - joint
        dir_joint /= np.linalg.norm(dir_joint)
        joint = mule - L * dir_joint
        trailer_dir = 0.85 * trailer_dir + 0.15 * dir_joint
        trailer_dir /= np.linalg.norm(trailer_dir)
        trailer = joint - L * trailer_dir
        mule_traj.append(mule.copy()); joint_traj.append(joint.copy()); trailer_traj.append(trailer.copy())
    return np.array(mule_traj), np.array(joint_traj), np.array(trailer_traj)

def simulate_three_point_body(path, L=1.5, joint_flex=0.3):
    mule_traj, joint_traj, trailer_traj = [], [], []
    trailer_dir = np.array([-1.0, 0.0])
    for i, mule in enumerate(path):
        prev = path[i-1] if i > 0 else mule
        mule_dir = mule - prev
        mule_dir = mule_dir / np.linalg.norm(mule_dir) if np.linalg.norm(mule_dir) > 1e-6 else np.array([1.0, 0.0])
        angle_offset = joint_flex * np.sin(i * 0.05)
        rot = np.array([[np.cos(angle_offset), -np.sin(angle_offset)], [np.sin(angle_offset), np.cos(angle_offset)]])
        joint_dir = rot @ (-mule_dir)
        joint = mule + L * joint_dir
        trailer_dir = 0.9 * trailer_dir + 0.1 * joint_dir
        trailer_dir /= np.linalg.norm(trailer_dir)
        trailer = joint + L * trailer_dir
        mule_traj.append(mule); joint_traj.append(joint); trailer_traj.append(trailer)
    return np.array(mule_traj), np.array(joint_traj), np.array(trailer_traj)

def compute_error(trailer_traj, path):
    return np.array([np.min(np.linalg.norm(path - p, axis=1)) for p in trailer_traj])

def animate_three_points(path, mule_traj, joint_traj, trailer_traj, title):
    fig, ax = plt.subplots()
    ax.plot(path[:, 0], path[:, 1], 'b--', label='Figure-8 Path')
    mule_point, = ax.plot([], [], 'ro', label='Mule')
    joint_point, = ax.plot([], [], 'yo', label='Joint')
    trailer_point, = ax.plot([], [], 'ko', label='Trailer')
    ax.legend(); ax.axis('equal'); ax.set_title(title)
    def init():
        mule_point.set_data([], []); joint_point.set_data([], []); trailer_point.set_data([], [])
        return mule_point, joint_point, trailer_point
    def update(i):
        mule_point.set_data([mule_traj[i, 0]], [mule_traj[i, 1]])
        joint_point.set_data([joint_traj[i, 0]], [joint_traj[i, 1]])
        trailer_point.set_data([trailer_traj[i, 0]], [trailer_traj[i, 1]])
        return mule_point, joint_point, trailer_point
    ani = animation.FuncAnimation(fig, update, frames=len(mule_traj), init_func=init, interval=30, blit=True)
    plt.show()

def figure_eight(t):
    return np.array([10 * np.sin(t), 10 * np.sin(t) * np.cos(t)])

path = generate_frenet_frame(figure_eight)
mule_nc, joint_nc, trailer_nc = simulate_three_point_body(path)
mule_mpc, joint_mpc, trailer_mpc = mpc_controller(path)
animate_three_points(path, mule_nc, joint_nc, trailer_nc, "Without MPC Controller")
animate_three_points(path, mule_mpc, joint_mpc, trailer_mpc, "With MPC Controller")
err_nc, err_mpc = compute_error(trailer_nc, path), compute_error(trailer_mpc, path)
plt.figure(); plt.plot(err_nc, 'r', label='No Controller'); plt.plot(err_mpc, 'g', label='MPC Controller')
plt.xlabel('Step'); plt.ylabel('Trailer Error'); plt.legend(); plt.title('Trailer Error Comparison'); plt.show()
