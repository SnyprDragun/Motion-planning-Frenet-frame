#!/Users/subhodeep/venv/bin/python

# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib import animation
# from scipy.optimize import minimize

# # Parameters
# L = 1.5
# D = 2.0
# v = 1.0
# dt = 0.1
# N = 10  # MPC horizon

# # Reference path: circle of radius 10
# def ref_path(t):
#     R = 10.0
#     omega = v / R
#     x = R * np.cos(omega * t)
#     y = R * np.sin(omega * t)
#     theta = np.arctan2(y, x) + np.pi / 2
#     return x, y, theta

# # Dynamics
# def step(x, u):
#     x_, y_, theta_, phi_ = x
#     omega = u

#     dx = v * np.cos(theta_)
#     dy = v * np.sin(theta_)
#     dtheta = omega
#     v_t = v * np.cos(theta_ - phi_) + L * omega * np.sin(theta_ - phi_)
#     dphi = (v * np.sin(theta_ - phi_) - L * omega * np.cos(theta_ - phi_)) / D

#     return np.array([x_ + dx * dt, y_ + dy * dt, theta_ + dtheta * dt, phi_ + dphi * dt])

# # MPC cost function
# def mpc_cost(omega_seq, state0, t0):
#     state = np.array(state0)
#     cost = 0
#     for i in range(N):
#         omega = omega_seq[i]
#         state = step(state, omega)
#         x_ref, y_ref, _ = ref_path(t0 + (i+1)*dt)
#         cost += (state[0] - x_ref)**2 + (state[1] - y_ref)**2
#     return cost

# # Control using MPC
# def get_mpc_control(state, t0):
#     init_guess = np.zeros(N)
#     bounds = [(-1.5, 1.5)] * N  # omega bounds
#     res = minimize(mpc_cost, init_guess, args=(state, t0), bounds=bounds, method='SLSQP')
#     return res.x[0] if res.success else 0.0

# # Compute hitch and trailer positions
# def compute_hitch_trailer(state):
#     x, y, theta, phi = state
#     x_h = x - L * np.cos(theta)
#     y_h = y - L * np.sin(theta)
#     x_t = x_h - D * np.cos(phi)
#     y_t = y_h - D * np.sin(phi)
#     return (x_h, y_h), (x_t, y_t)

# # ---- Simulation ----
# T = 100
# steps = int(T / dt)
# states = []
# hitches = []
# trailers = []
# state = [10, 0, np.pi/2, np.pi/2]  # Start on circle
# times = []

# for i in range(steps):
#     t = i * dt
#     omega = get_mpc_control(state, t)
#     state = step(state, omega)
#     states.append(state)
#     hitch, trailer = compute_hitch_trailer(state)
#     hitches.append(hitch)
#     trailers.append(trailer)
#     times.append(t)

# states = np.array(states)
# hitches = np.array(hitches)
# trailers = np.array(trailers)

# # ---- Animation ----
# fig, ax = plt.subplots()
# ax.set_aspect('equal')
# ax.set_xlim(-12, 12)
# ax.set_ylim(-12, 12)

# line_path, = ax.plot([], [], 'b--', label='Reference Path')
# line_traj, = ax.plot([], [], 'r-', label='Mule Path')
# mule_point, = ax.plot([], [], 'ro', label='Mule')
# hitch_point, = ax.plot([], [], 'yo', label='Hitch')
# trailer_point, = ax.plot([], [], 'ko', label='Trailer')
# link1, = ax.plot([], [], 'r-', lw=1.5)
# link2, = ax.plot([], [], 'k-', lw=1.5)

# circle = plt.Circle((0, 0), 10, color='gray', fill=False, linestyle='--')
# ax.add_patch(circle)
# ax.legend()

# def init():
#     line_path.set_data([], [])
#     line_traj.set_data([], [])
#     mule_point.set_data([], [])
#     hitch_point.set_data([], [])
#     trailer_point.set_data([], [])
#     link1.set_data([], [])
#     link2.set_data([], [])
#     return line_path, line_traj, mule_point, hitch_point, trailer_point, link1, link2

# def update(i):
#     mule = states[i, :2]
#     hitch = hitches[i]
#     trailer = trailers[i]

#     # Mule trajectory so far
#     line_traj.set_data(states[:i+1, 0], states[:i+1, 1])
#     mule_point.set_data([mule[0]], [mule[1]])
#     hitch_point.set_data([hitch[0]], [hitch[1]])
#     trailer_point.set_data([trailer[0]], [trailer[1]])
#     link1.set_data([mule[0], hitch[0]], [mule[1], hitch[1]])
#     link2.set_data([hitch[0], trailer[0]], [hitch[1], trailer[1]])
#     return line_traj, mule_point, hitch_point, trailer_point, link1, link2

# ani = animation.FuncAnimation(fig, update, frames=steps, init_func=init, interval=50, blit=True)
# plt.title("MPC Tracking a Circle with Hitch & Trailer")
# plt.show()



import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.optimize import minimize
import math

# ===================================
#  Class: CartesianFrenetConverter
# ===================================
class CartesianFrenetConverter:
    """Convert states between Cartesian and Frenet coordinate systems"""

    @staticmethod
    def cartesian_to_frenet(rs, rx, ry, rtheta, rkappa, rdkappa, x, y, v, a, theta, kappa):
        dx = x - rx
        dy = y - ry
        cos_r = math.cos(rtheta)
        sin_r = math.sin(rtheta)
        cross = cos_r * dy - sin_r * dx
        d = math.copysign(math.hypot(dx, dy), cross)
        delta_theta = theta - rtheta
        tan_delta = math.tan(delta_theta)
        cos_delta = math.cos(delta_theta)
        one_minus_kr_d = 1 - rkappa * d
        d_dot = one_minus_kr_d * tan_delta
        kappa_r_d_prime = rdkappa * d + rkappa * d_dot
        d_ddot = (-kappa_r_d_prime * tan_delta +
                  one_minus_kr_d / (cos_delta * cos_delta) *
                  (kappa * one_minus_kr_d / cos_delta - rkappa))
        s = rs
        s_dot = v * cos_delta / one_minus_kr_d
        delta_theta_prime = one_minus_kr_d / cos_delta * kappa - rkappa
        s_ddot = (a * cos_delta - s_dot * s_dot *
                  (d_dot * delta_theta_prime - kappa_r_d_prime)) / one_minus_kr_d
        return [s, s_dot, s_ddot], [d, d_dot, d_ddot]

# ===========================
#  Class: ParametricFunctions
# ===========================
class ParametricFunctions:
    @staticmethod
    def circle(t, R=10.0, v=1.0):
        omega = v / R
        x = R * np.cos(omega * t)
        y = R * np.sin(omega * t)
        theta = np.arctan2(y, x) + np.pi / 2
        return x, y, theta

    @staticmethod
    def ellipse(t, a=12.0, b=8.0, v=1.0):
        omega = v / ((a + b) / 2.0)
        x = a * np.cos(omega * t)
        y = b * np.sin(omega * t)
        theta = np.arctan2(b * np.cos(omega * t), -a * np.sin(omega * t))
        return x, y, theta

# =====================
#  Class: Controller
# =====================
class Controller:
    def __init__(self, N=10, dt=0.1, v=1.0):
        self.N = N
        self.dt = dt
        self.v = v

    def step_dynamics(self, state, omega, L, D):
        x, y, theta, phi = state
        dx = self.v * np.cos(theta)
        dy = self.v * np.sin(theta)
        dtheta = omega
        dphi = (self.v * np.sin(theta - phi) - L * omega * np.cos(theta - phi)) / D
        return np.array([x + dx * self.dt, y + dy * self.dt, theta + dtheta * self.dt, phi + dphi * self.dt])

    def cost(self, omega_seq, state0, t0, path_func, L, D):
        state = np.array(state0)
        J = 0
        for i in range(self.N):
            state = self.step_dynamics(state, omega_seq[i], L, D)
            x_ref, y_ref, _ = path_func(t0 + (i + 1) * self.dt)
            J += (state[0] - x_ref) ** 2 + (state[1] - y_ref) ** 2
        return J

    def mpc(self, x, y, theta, phi, path_func, t0=0, L=1.5, D=2.0):
        state0 = [x, y, theta, phi]
        init_guess = np.zeros(self.N)
        bounds = [(-1.5, 1.5)] * self.N
        res = minimize(self.cost, init_guess, args=(state0, t0, path_func, L, D),
                       bounds=bounds, method='SLSQP')
        omega_opt = res.x[0] if res.success else 0.0
        return omega_opt

# =====================
#  Class: PathPlanningFF
# =====================
class PathPlanningFF:
    def __init__(self, path_func, controller, L=1.5, D=2.0, dt=0.1, T=100):
        self.path_func = path_func
        self.controller = controller
        self.L = L
        self.D = D
        self.dt = dt
        self.T = T
        self.steps = int(T / dt)

        self.x, self.y, self.theta = 10, 0, np.pi / 2
        self.phi = np.pi / 2
        self.v = controller.v

        self.states = []
        self.hitches = []
        self.trailers = []
        self.frenet_states = []
        self.frenet_hitches = []
        self.frenet_trailers = []

        # Precompute Frenet reference path
        self.ref_s = []
        self.ref_d = []

    def compute_hitch_trailer(self):
        x_h = self.x - self.L * np.cos(self.theta)
        y_h = self.y - self.L * np.sin(self.theta)
        x_t = x_h - self.D * np.cos(self.phi)
        y_t = y_h - self.D * np.sin(self.phi)
        return np.array([x_h, y_h]), np.array([x_t, y_t])

    def simulate(self):
        for i in range(self.steps):
            t = i * self.dt
            omega = self.controller.mpc(self.x, self.y, self.theta, self.phi,
                                        self.path_func, t0=t, L=self.L, D=self.D)

            # dynamics
            self.phi += ((self.v * np.sin(self.theta - self.phi) - self.L * omega * np.cos(self.theta - self.phi)) / self.D) * self.dt
            self.x += self.v * np.cos(self.theta) * self.dt
            self.y += self.v * np.sin(self.theta) * self.dt
            self.theta += omega * self.dt

            hitch, trailer = self.compute_hitch_trailer()

            # Log Cartesian
            self.states.append([self.x, self.y])
            self.hitches.append(hitch)
            self.trailers.append(trailer)

            # Log Frenet (for mule, hitch, trailer)
            x_ref, y_ref, theta_ref = self.path_func(t)
            for point, storage in zip([[self.x, self.y, self.theta],
                                       [hitch[0], hitch[1], self.theta],
                                       [trailer[0], trailer[1], self.theta]],
                                      [self.frenet_states, self.frenet_hitches, self.frenet_trailers]):
                s_cond, d_cond = CartesianFrenetConverter.cartesian_to_frenet(
                    rs=t*self.v, rx=x_ref, ry=y_ref, rtheta=theta_ref,
                    rkappa=0.0, rdkappa=0.0, x=point[0], y=point[1],
                    v=self.v, a=0.0, theta=point[2], kappa=0.0
                )
                storage.append([s_cond[0], d_cond[0]])

            # store Frenet reference path
            self.ref_s.append(t * self.v)
            self.ref_d.append(0.0)

        # convert to arrays
        self.states = np.array(self.states)
        self.hitches = np.array(self.hitches)
        self.trailers = np.array(self.trailers)
        self.frenet_states = np.array(self.frenet_states)
        self.frenet_hitches = np.array(self.frenet_hitches)
        self.frenet_trailers = np.array(self.frenet_trailers)
        self.ref_s = np.array(self.ref_s)
        self.ref_d = np.array(self.ref_d)

    def animate(self):
        fig, (ax_cart, ax_frenet) = plt.subplots(1, 2, figsize=(12, 6))

        # ==== Cartesian Plot ====
        ax_cart.set_aspect('equal')
        ax_cart.set_xlim(-12, 12)
        ax_cart.set_ylim(-12, 12)
        circle = plt.Circle((0, 0), 10, color='gray', fill=False, linestyle='--')
        ax_cart.add_patch(circle)

        line_traj, = ax_cart.plot([], [], 'r-', label='Mule Path')
        mule_point, = ax_cart.plot([], [], 'ro', label='Mule')
        hitch_point, = ax_cart.plot([], [], 'yo', label='Hitch')
        trailer_point, = ax_cart.plot([], [], 'ko', label='Trailer')
        link1, = ax_cart.plot([], [], 'r-', lw=1.5)
        link2, = ax_cart.plot([], [], 'k-', lw=1.5)
        ax_cart.legend()
        ax_cart.set_title("Cartesian Frame")

        # ==== Frenet Plot ====
        ax_frenet.set_xlim(0, self.steps * self.dt * self.v)
        ax_frenet.set_ylim(-5, 5)
        ax_frenet.axhline(0, color='gray', linestyle='--', label='Reference Path (d=0)')

        line_frenet, = ax_frenet.plot([], [], 'r-', label='Mule s-d')
        mule_frenet_point, = ax_frenet.plot([], [], 'ro', label='Mule')
        ax_frenet.legend()
        ax_frenet.set_xlabel("s (arc length)")
        ax_frenet.set_ylabel("d (lateral deviation)")
        ax_frenet.set_title("Frenet Frame")

        # ==== Precompute Frenet Coordinates ====
        frenet_coords = []
        s_accum = 0
        for i in range(len(self.states)):
            x, y = self.states[i]
            x_ref, y_ref, theta_ref = self.path_func(i * self.dt)
            # project mule on ref path: s increases linearly with v*dt
            s_accum += self.v * self.dt
            d = np.sin(theta_ref) * (x - x_ref) - np.cos(theta_ref) * (y - y_ref)
            frenet_coords.append([s_accum, d])
        frenet_coords = np.array(frenet_coords)

        def init():
            # Cartesian
            line_traj.set_data([], [])
            mule_point.set_data([], [])
            hitch_point.set_data([], [])
            trailer_point.set_data([], [])
            link1.set_data([], [])
            link2.set_data([], [])
            # Frenet
            line_frenet.set_data([], [])
            mule_frenet_point.set_data([], [])
            return (line_traj, mule_point, hitch_point, trailer_point, link1, link2,
                    line_frenet, mule_frenet_point)

        def update(i):
            # === Cartesian ===
            mule = self.states[i]
            hitch = self.hitches[i]
            trailer = self.trailers[i]

            line_traj.set_data(self.states[:i+1, 0], self.states[:i+1, 1])
            mule_point.set_data([mule[0]], [mule[1]])
            hitch_point.set_data([hitch[0]], [hitch[1]])
            trailer_point.set_data([trailer[0]], [trailer[1]])
            link1.set_data([mule[0], hitch[0]], [mule[1], hitch[1]])
            link2.set_data([hitch[0], trailer[0]], [hitch[1], trailer[1]])

            # === Frenet === (Only Mule)
            line_frenet.set_data(frenet_coords[:i+1, 0], frenet_coords[:i+1, 1])
            mule_frenet_point.set_data([frenet_coords[i, 0]], [frenet_coords[i, 1]])

            return (line_traj, mule_point, hitch_point, trailer_point, link1, link2,
                    line_frenet, mule_frenet_point)

        ani = animation.FuncAnimation(fig, update, frames=self.steps,
                                    init_func=init, interval=50, blit=True)
        plt.show()

# ============
#   MAIN
# ============
if __name__ == "__main__":
    path_func = lambda t: ParametricFunctions.circle(t, R=10.0, v=1.0)
    controller = Controller(N=10, dt=0.1, v=1.0)
    sim = PathPlanningFF(path_func, controller)
    sim.simulate()
    sim.animate()
