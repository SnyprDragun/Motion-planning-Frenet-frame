#!/Users/subhodeep/venv/bin/python
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.optimize import minimize

class CartesianFrenetConverter:
    """
    A class for converting states between Cartesian and Frenet coordinate systems
    """

    @ staticmethod
    def cartesian_to_frenet(rs, rx, ry, rtheta, rkappa, rdkappa, x, y, v, a, theta, kappa):
        """
        Convert state from Cartesian coordinate to Frenet coordinate

        Parameters
        ----------
            rs: reference line s-coordinate
            rx, ry: reference point coordinates
            rtheta: reference point heading
            rkappa: reference point curvature
            rdkappa: reference point curvature rate
            x, y: current position
            v: velocity
            a: acceleration
            theta: heading angle
            kappa: curvature

        Returns
        -------
            s_condition: [s(t), s'(t), s''(t)]
            d_condition: [d(s), d'(s), d''(s)]
        """
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

    @ staticmethod
    def frenet_to_cartesian(rs, rx, ry, rtheta, rkappa, rdkappa, s_condition, d_condition):
        """
        Convert state from Frenet coordinate to Cartesian coordinate

        Parameters
        ----------
            rs: reference line s-coordinate
            rx, ry: reference point coordinates
            rtheta: reference point heading
            rkappa: reference point curvature
            rdkappa: reference point curvature rate
            s_condition: [s(t), s'(t), s''(t)]
            d_condition: [d(s), d'(s), d''(s)]

        Returns
        -------
            x, y: position
            theta: heading angle
            kappa: curvature
            v: velocity
            a: acceleration
        """
        if abs(rs - s_condition[0]) >= 1.0e-6:
            raise ValueError(
                "The reference point s and s_condition[0] don't match")

        cos_theta_r = math.cos(rtheta)
        sin_theta_r = math.sin(rtheta)

        x = rx - sin_theta_r * d_condition[0]
        y = ry + cos_theta_r * d_condition[0]

        one_minus_kappa_r_d = 1 - rkappa * d_condition[0]

        tan_delta_theta = d_condition[1] / one_minus_kappa_r_d
        delta_theta = math.atan2(d_condition[1], one_minus_kappa_r_d)
        cos_delta_theta = math.cos(delta_theta)

        theta = CartesianFrenetConverter.normalize_angle(delta_theta + rtheta)

        kappa_r_d_prime = rdkappa * d_condition[0] + rkappa * d_condition[1]

        kappa = (((d_condition[2] + kappa_r_d_prime * tan_delta_theta) *
                  cos_delta_theta * cos_delta_theta) / one_minus_kappa_r_d + rkappa) * \
            cos_delta_theta / one_minus_kappa_r_d

        d_dot = d_condition[1] * s_condition[1]
        v = math.sqrt(one_minus_kappa_r_d * one_minus_kappa_r_d *
                      s_condition[1] * s_condition[1] + d_dot * d_dot)

        delta_theta_prime = one_minus_kappa_r_d / cos_delta_theta * kappa - rkappa

        a = (s_condition[2] * one_minus_kappa_r_d / cos_delta_theta +
             s_condition[1] * s_condition[1] / cos_delta_theta *
             (d_condition[1] * delta_theta_prime - kappa_r_d_prime))

        return x, y, theta, kappa, v, a

    @ staticmethod
    def normalize_angle(angle):
        """
        Normalize angle to [-pi, pi]
        """
        a = math.fmod(angle + math.pi, 2.0 * math.pi)
        if a < 0.0:
            a += 2.0 * math.pi
        return a - math.pi


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

    @staticmethod
    def figure_eight(t, a=10.0, v=1.0):
        omega = v / a
        x = a * np.sin(omega * t)
        y = a * np.sin(omega * t) * np.cos(omega * t)
        dx = a * omega * np.cos(omega * t)
        dy = a * omega * np.cos(2 * omega * t)
        theta = np.arctan2(dy, dx)
        return x, y, theta


class Controller:
    def __init__(self, N=12, dt=0.1, v=1.0,
                 w_d_mule=1200.0, w_d_hitch=0.0, w_d_trailer=0.0,
                 w_omega=0.01, w_omega_rate=1.0,
                 omega_bounds=(-2.0, 2.0), max_delta_omega=0.5,
                 use_exp_smooth=False, smooth_alpha=0.3):
        """
        MPC controller using only lateral Frenet errors of mule, hitch, trailer.
        - N: horizon steps
        - dt: timestep
        - v: forward speed (constant)
        - w_d_mule/hitch/trailer: weights on lateral deviations
        - w_omega, w_omega_rate: regularization on omega magnitude and rate
        - omega_bounds: bounds on steering rate omega
        - max_delta_omega: max allowed change in omega per step
        - use_exp_smooth: apply exponential smoothing on omega output
        """
        self.N = N
        self.dt = dt
        self.v = v

        self.w_d_mule = w_d_mule
        self.w_d_hitch = w_d_hitch
        self.w_d_trailer = w_d_trailer

        self.w_omega = w_omega
        self.w_omega_rate = w_omega_rate

        self.omega_min, self.omega_max = omega_bounds
        self.prev_omega = 0.0
        self.max_delta_omega = max_delta_omega

        self.use_exp_smooth = use_exp_smooth
        self.smooth_alpha = smooth_alpha

    def wrap_angle(self, angle):
        """Wrap angle to [-pi, pi]."""
        return (angle + np.pi) % (2 * np.pi) - np.pi

    def step_dynamics(self, state, omega, L, D):
        x, y, theta, phi = state
        dx = self.v * np.cos(theta)
        dy = self.v * np.sin(theta)
        dtheta = omega
        dphi = (self.v * np.sin(theta - phi) - L * omega * np.cos(theta - phi)) / D

        theta_new = self.wrap_angle(theta + dtheta * self.dt)
        phi_new = self.wrap_angle(phi + dphi * self.dt)

        return np.array([
            x + dx * self.dt,
            y + dy * self.dt,
            theta_new,
            phi_new
        ])

    def lateral_deviation(self, x, y, theta_ref, rx, ry):
        """
        Compute lateral deviation d for pose (x,y) relative to reference path point (rx, ry, theta_ref).
        """
        return math.sin(theta_ref) * (x - rx) - math.cos(theta_ref) * (y - ry)

    def compute_hitch_trailer_pose(self, mule_pose, phi, L, D):
        """
        Compute hitch and trailer poses from mule pose and hitch angle phi.
        Returns: (x_hitch, y_hitch, theta_hitch), (x_trailer, y_trailer, theta_trailer)
        """
        x, y, theta = mule_pose

        # Hitch position behind mule by L
        x_hitch = x - L * np.cos(theta)
        y_hitch = y - L * np.sin(theta)
        theta_hitch = self.wrap_angle(theta + phi)

        # Trailer position behind hitch by D
        x_trailer = x_hitch - D * np.cos(theta_hitch)
        y_trailer = y_hitch - D * np.sin(theta_hitch)
        theta_trailer = theta_hitch  # trailer heading same as hitch

        return (x_hitch, y_hitch, theta_hitch), (x_trailer, y_trailer, theta_trailer)

    def cost(self, omega_seq, state0, t0, path_func, L, D):
        state = np.array(state0, dtype=float)
        J = 0.0
        prev_omega = self.prev_omega

        for i in range(self.N):
            omega = float(omega_seq[i])
            state = self.step_dynamics(state, omega, L, D)
            x, y, theta, phi = state

            t_pred = t0 + (i + 1) * self.dt

            # Reference point for mule at current predicted time
            rx, ry, theta_ref = path_func(t_pred)

            # Time-shifted reference points for hitch and trailer (accounting for physical offset)
            t_hitch = max(t_pred - L / self.v, 0.0)
            t_trailer = max(t_pred - (L + D) / self.v, 0.0)

            rx_hitch, ry_hitch, theta_ref_hitch = path_func(t_hitch)
            rx_trailer, ry_trailer, theta_ref_trailer = path_func(t_trailer)

            # Lateral deviation for mule
            d_mule = self.lateral_deviation(x, y, theta_ref, rx, ry)

            # Hitch and trailer poses
            (x_h, y_h, th_h), (x_t, y_t, th_t) = self.compute_hitch_trailer_pose((x, y, theta), phi, L, D)

            # Lateral deviations for hitch and trailer relative to their shifted reference points
            d_hitch = self.lateral_deviation(x_h, y_h, theta_ref_hitch, rx_hitch, ry_hitch)
            d_trailer = self.lateral_deviation(x_t, y_t, theta_ref_trailer, rx_trailer, ry_trailer)

            # Accumulate weighted lateral deviation costs
            J += self.w_d_mule * (d_mule ** 2)
            J += self.w_d_hitch * (d_hitch ** 2)
            J += self.w_d_trailer * (d_trailer ** 2)

            # Regularize omega and omega rate change
            J += self.w_omega * (omega ** 2)
            J += self.w_omega_rate * ((omega - prev_omega) ** 2)

            heading_err = math.atan2(math.sin(theta - theta_ref), math.cos(theta - theta_ref))
            J += 200.0 * (heading_err ** 2)  # small weight, tune as needed

            prev_omega = omega

        return J

    def mpc(self, x, y, theta, phi, path_func, t0=0.0, L=1.5, D=2.0):
        state0 = [x, y, theta, phi]
        init_guess = np.full(self.N, self.prev_omega, dtype=float)
        bounds = [(self.omega_min, self.omega_max)] * self.N

        res = minimize(self.cost, init_guess,
                       args=(state0, t0, path_func, L, D),
                       bounds=bounds, method='SLSQP',
                       options={'maxiter': 200, 'ftol': 1e-4})

        if res.success:
            omega_raw = float(res.x[0])
        else:
            omega_raw = self.prev_omega

        if self.use_exp_smooth:
            omega_out = self.smooth_alpha * omega_raw + (1.0 - self.smooth_alpha) * self.prev_omega
        else:
            omega_out = omega_raw

        maxd = self.max_delta_omega
        omega_out = float(np.clip(omega_out, self.prev_omega - maxd, self.prev_omega + maxd))
        self.prev_omega = omega_out
        return omega_out


class PathPlanningFrenetFrame:
    def __init__(self, path_func, start, controller, L=1.5, D=2.0, dt=0.1, T=60):
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


if __name__ == "__main__":
    path_func = lambda t: ParametricFunctions.figure_eight(t, a=10.0, v=1.0)
    controller = Controller(N=10, dt=0.1, v=1.0)
    sim = PathPlanningFrenetFrame(path_func, (0, 0, np.pi/4), controller)
    sim.simulate()
    sim.animate()
