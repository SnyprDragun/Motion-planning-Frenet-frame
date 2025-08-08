#!/Users/subhodeep/venv/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from CartesianFrenetConverter import CartesianFrenetConverter

class PathPlanningFrenetFrame:
    def __init__(self, path_func, start, controller, L=1.5, D=2.0, dt=0.1, T=100):
        self.path_func = path_func
        self.controller = controller
        self.L = L
        self.D = D
        self.dt = dt
        self.T = T
        self.steps = int(T / dt)
        self.start = start
        self.x, self.y, self.theta = self.start
        self.phi = self.theta
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
