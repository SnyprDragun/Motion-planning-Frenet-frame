#!/Users/subhodeep/venv/bin/python
import numpy as np
import matplotlib.pyplot as plt

# # ================== Helper Functions ==================
# def generate_figure_eight(num_points=500, radius=10):
#     t = np.linspace(0, 2*np.pi, num_points)
#     x = radius * np.sin(t)
#     y = radius * np.sin(t) * np.cos(t)
#     return np.vstack((x, y)).T

# def cartesian_to_frenet(path):
#     s = [0]
#     for i in range(1, len(path)):
#         ds = np.linalg.norm(path[i] - path[i-1])
#         s.append(s[-1] + ds)
#     s = np.array(s)
#     d = np.zeros_like(s)
#     return np.vstack((s, d)).T

# def plot_path(cartesian, frenet=None, title="Path"):
#     plt.figure()
#     plt.plot(cartesian[:,0], cartesian[:,1], label='Cartesian')
#     if frenet is not None:
#         plt.plot(frenet[:,0], frenet[:,1], label='Frenet')
#     plt.legend()
#     plt.axis('equal')
#     plt.title(title)
#     plt.show()

# # ================== Controllers ==================
# def pure_pursuit(current_pos, path, lookahead=1.0):
#     dists = np.linalg.norm(path - current_pos, axis=1)
#     idx = np.argmin(dists)
#     goal_idx = min(idx + 5, len(path)-1)
#     goal = path[goal_idx]
#     direction = (goal - current_pos) / np.linalg.norm(goal - current_pos)
#     next_pos = current_pos + 0.2 * direction
#     return next_pos

# # ================== Simulation Modules ==================
# def simulate_tracking(path):
#     pos = path[0].copy()
#     traj = [pos]
#     for _ in range(1000):
#         pos = pure_pursuit(pos, path)
#         traj.append(pos)
#         if np.linalg.norm(pos - path[-1]) < 0.5:
#             break
#     return np.array(traj)

# # ================== Part 1: Figure of Eight ==================
# path_cart = generate_figure_eight()
# path_frenet = cartesian_to_frenet(path_cart)
# plot_path(path_cart, path_frenet, "Figure-8 in Cartesian and Frenet")

# # ================== Part 2: Tracking ==================
# traj = simulate_tracking(path_cart)
# plot_path(path_cart, traj, "Tracking Figure-8")

# # ================== Part 3: Obstacle Avoidance ==================
# obstacle = np.array([0, 0])
# path_avoid = np.array([p + [0, 3] if np.linalg.norm(p - obstacle) < 5 else p for p in path_cart])
# traj_avoid = simulate_tracking(path_avoid)
# plot_path(path_avoid, traj_avoid, "Obstacle Avoidance")

# # ================== Part 3b and 4: Trailer Kinematics ==================
# def trailer_kinematics(pos, theta, phi, v=0.2, L=2.0, D=2.0, dt=0.1):
#     x, y = pos
#     x += v * np.cos(theta) * dt
#     y += v * np.sin(theta) * dt
#     theta += (v / L) * np.tan(phi) * dt
#     return np.array([x, y]), theta

# def simulate_trailer(path):
#     pos = path[0].copy()
#     theta = 0.0
#     phi = 0.0
#     traj = [pos]
#     for _ in range(1000):
#         pos, theta = trailer_kinematics(pos, theta, phi)
#         traj.append(pos)
#         if np.linalg.norm(pos - path[-1]) < 0.5:
#             break
#     return np.array(traj)

# traj_trailer = simulate_trailer(path_cart)
# plot_path(path_cart, traj_trailer, "Trailer Tracking Figure-8")

# # ================== Part 4 Reverse Tracking ==================
# traj_reverse = np.flip(traj_trailer, axis=0)
# plot_path(path_cart, traj_reverse, "Reverse Trailer Tracking")

# print("Simulation complete.")



import numpy as np
import matplotlib.pyplot as plt
import math

# ================== Frenet Conversion Function ==================
def cartesian_to_frenet(rs, rx, ry, rtheta, rkappa, rdkappa, x, y, v, a, theta, kappa):
    dx = x - rx
    dy = y - ry

    cos_theta_r = math.cos(rtheta)
    sin_theta_r = math.sin(rtheta)

    cross_rd_nd = cos_theta_r * dy - sin_theta_r * dx
    d = math.copysign(math.hypot(dx, dy), cross_rd_nd)

    delta_theta = theta - rtheta
    tan_delta_theta = math.tan(delta_theta)
    cos_delta_theta = math.cos(delta_theta)

    one_minus_kappa_r_d = 1 - rkappa * d
    d_dot = one_minus_kappa_r_d * tan_delta_theta

    kappa_r_d_prime = rdkappa * d + rkappa * d_dot

    d_ddot = (-kappa_r_d_prime * tan_delta_theta +
              one_minus_kappa_r_d / (cos_delta_theta * cos_delta_theta) *
              (kappa * one_minus_kappa_r_d / cos_delta_theta - rkappa))

    s = rs
    s_dot = v * cos_delta_theta / one_minus_kappa_r_d

    delta_theta_prime = one_minus_kappa_r_d / cos_delta_theta * kappa - rkappa
    s_ddot = (a * cos_delta_theta -
              s_dot * s_dot *
              (d_dot * delta_theta_prime - kappa_r_d_prime)) / one_minus_kappa_r_d

    return [s, s_dot, s_ddot], [d, d_dot, d_ddot]

# ================== Parametric Curve (Figure Eight) ==================
def figure_eight(t):
    x = 10 * np.sin(t)
    y = 10 * np.sin(t) * np.cos(t)
    return x, y

# ================== Compute Frenet Frame ==================
def compute_frenet_for_curve(num_points=300):
    t_values = np.linspace(0, 2*np.pi, num_points)
    points = np.array([figure_eight(t) for t in t_values])

    s_vals = [0.0]
    d_vals = []
    for i in range(1, len(points)):
        s_vals.append(s_vals[-1] + np.linalg.norm(points[i] - points[i-1]))

    s_vals = np.array(s_vals)

    # Approximate tangent angle rtheta for the path
    thetas = np.arctan2(np.gradient(points[:, 1]), np.gradient(points[:, 0]))

    # For simplicity, assume rkappa and rdkappa are zero (straight line segments)
    frenet_s = []
    frenet_d = []

    for i in range(len(points)):
        rs = s_vals[i]
        rx, ry = points[i]
        rtheta = thetas[i]
        rkappa = 0.0
        rdkappa = 0.0

        # For the test, we use the same point as vehicle state with same heading
        s_cond, d_cond = cartesian_to_frenet(rs, rx, ry, rtheta, rkappa, rdkappa,
                                             rx, ry, 1.0, 0.0, rtheta, 0.0)
        frenet_s.append(s_cond[0])
        frenet_d.append(d_cond[0])

    return points[:, 0], points[:, 1], np.array(frenet_s), np.array(frenet_d)

# ================== Main Execution ==================
x, y, s, d = compute_frenet_for_curve()

# ================== Plotting ==================
fig, ax = plt.subplots(1, 2, figsize=(12, 5))

# Cartesian Frame
ax[0].plot(x, y, 'b')
ax[0].set_title("Curve in Cartesian Frame (x vs y)")
ax[0].set_xlabel("x")
ax[0].set_ylabel("y")
ax[0].axis('equal')

# Frenet Frame
ax[1].plot(s, d, 'r')
ax[1].set_title("Curve in Frenet Frame (s vs d)")
ax[1].set_xlabel("s (arc length)")
ax[1].set_ylabel("d (lateral offset)")

plt.tight_layout()
plt.show()
