#!/Users/subhodeep/venv/bin/python
import numpy as np
import matplotlib.pyplot as plt

# ================== Helper Functions ==================
def generate_figure_eight(num_points=500, radius=10):
    t = np.linspace(0, 2*np.pi, num_points)
    x = radius * np.sin(t)
    y = radius * np.sin(t) * np.cos(t)
    return np.vstack((x, y)).T

def cartesian_to_frenet(path):
    s = [0]
    for i in range(1, len(path)):
        ds = np.linalg.norm(path[i] - path[i-1])
        s.append(s[-1] + ds)
    s = np.array(s)
    d = np.zeros_like(s)
    return np.vstack((s, d)).T

def plot_path(cartesian, frenet=None, title="Path"):
    plt.figure()
    plt.plot(cartesian[:,0], cartesian[:,1], label='Cartesian')
    if frenet is not None:
        plt.plot(frenet[:,0], frenet[:,1], label='Frenet')
    plt.legend()
    plt.axis('equal')
    plt.title(title)
    plt.show()

# ================== Controllers ==================
def pure_pursuit(current_pos, path, lookahead=1.0):
    dists = np.linalg.norm(path - current_pos, axis=1)
    idx = np.argmin(dists)
    goal_idx = min(idx + 5, len(path)-1)
    goal = path[goal_idx]
    direction = (goal - current_pos) / np.linalg.norm(goal - current_pos)
    next_pos = current_pos + 0.2 * direction
    return next_pos

# ================== Simulation Modules ==================
def simulate_tracking(path):
    pos = path[0].copy()
    traj = [pos]
    for _ in range(1000):
        pos = pure_pursuit(pos, path)
        traj.append(pos)
        if np.linalg.norm(pos - path[-1]) < 0.5:
            break
    return np.array(traj)

# ================== Part 1: Figure of Eight ==================
path_cart = generate_figure_eight()
path_frenet = cartesian_to_frenet(path_cart)
plot_path(path_cart, path_frenet, "Figure-8 in Cartesian and Frenet")

# ================== Part 2: Tracking ==================
traj = simulate_tracking(path_cart)
plot_path(path_cart, traj, "Tracking Figure-8")

# ================== Part 3: Obstacle Avoidance ==================
obstacle = np.array([0, 0])
path_avoid = np.array([p + [0, 3] if np.linalg.norm(p - obstacle) < 5 else p for p in path_cart])
traj_avoid = simulate_tracking(path_avoid)
plot_path(path_avoid, traj_avoid, "Obstacle Avoidance")

# ================== Part 3b and 4: Trailer Kinematics ==================
def trailer_kinematics(pos, theta, phi, v=0.2, L=2.0, D=2.0, dt=0.1):
    x, y = pos
    x += v * np.cos(theta) * dt
    y += v * np.sin(theta) * dt
    theta += (v / L) * np.tan(phi) * dt
    return np.array([x, y]), theta

def simulate_trailer(path):
    pos = path[0].copy()
    theta = 0.0
    phi = 0.0
    traj = [pos]
    for _ in range(1000):
        pos, theta = trailer_kinematics(pos, theta, phi)
        traj.append(pos)
        if np.linalg.norm(pos - path[-1]) < 0.5:
            break
    return np.array(traj)

traj_trailer = simulate_trailer(path_cart)
plot_path(path_cart, traj_trailer, "Trailer Tracking Figure-8")

# ================== Part 4 Reverse Tracking ==================
traj_reverse = np.flip(traj_trailer, axis=0)
plot_path(path_cart, traj_reverse, "Reverse Trailer Tracking")

print("Simulation complete.")
