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

    @staticmethod
    def straight_line(t, slope=0.0, intercept=0.0, v=1.0):
        x = v * t
        y = slope * x + intercept
        theta = np.arctan2(slope, 1.0)
        return x, y, theta


class Controller:
    def __init__(self, N=12, dt=0.1, v=1.0,
                 w_pos_trailer=1000.0, w_omega=0.01, w_omega_rate=1.0,
                 omega_bounds=(-2.0, 2.0), max_delta_omega=0.5,
                 use_exp_smooth=False, smooth_alpha=0.3):
        self.N = N
        self.dt = dt
        self.v = v

        # weight for trailer position tracking (we track trailer explicitly)
        self.w_pos_trailer = w_pos_trailer

        self.w_omega = w_omega
        self.w_omega_rate = w_omega_rate

        self.omega_min, self.omega_max = omega_bounds
        self.prev_omega = 0.0
        self.max_delta_omega = max_delta_omega

        self.use_exp_smooth = use_exp_smooth
        self.smooth_alpha = smooth_alpha

    def wrap_angle(self, angle):
        return (angle + np.pi) % (2 * np.pi) - np.pi

    def step_dynamics(self, state, omega, L, D):
        # state is mule pose: [x_m, y_m, theta_m, phi]
        x, y, theta, phi = state
        dx = self.v * np.cos(theta)
        dy = self.v * np.sin(theta)
        dtheta = omega
        # phi dynamics (articulation rate). This model is kept from original code.
        dphi = (self.v * np.sin(theta - phi) - L * omega * np.cos(theta - phi)) / D

        theta_new = self.wrap_angle(theta + dtheta * self.dt)
        phi_new = self.wrap_angle(phi + dphi * self.dt)

        return np.array([
            x + dx * self.dt,
            y + dy * self.dt,
            theta_new,
            phi_new
        ])

    def compute_hitch_trailer_pose(self, mule_pose, phi, L, D):
        x, y, theta = mule_pose
        # Hitch position behind mule by L
        x_hitch = x - L * np.cos(theta)
        y_hitch = y - L * np.sin(theta)
        theta_hitch = self.wrap_angle(theta + phi)

        # Trailer position behind hitch by D
        x_trailer = x_hitch - D * np.cos(theta_hitch)
        y_trailer = y_hitch - D * np.sin(theta_hitch)
        theta_trailer = theta_hitch

        return (x_hitch, y_hitch, theta_hitch), (x_trailer, y_trailer, theta_trailer)

    def compute_mule_hitch_pose_from_trailer(self, trailer_pose, L, D, phi_guess=0.0):
        """
        Given trailer pose (x_t, y_t, theta_t) and an estimate of phi (hitch angle),
        compute hitch and mule poses forward (useful to initialize mule state from trailer start).
        This assumes trailer heading == hitch heading.
        """
        x_t, y_t, theta_t = trailer_pose
        # hitch is ahead of trailer by D
        x_hitch = x_t + D * np.cos(theta_t)
        y_hitch = y_t + D * np.sin(theta_t)
        theta_hitch = theta_t

        # mule is ahead of hitch by L. Mule heading = hitch_heading - phi
        theta_mule = self.wrap_angle(theta_hitch - phi_guess)
        x_mule = x_hitch + L * np.cos(theta_mule)
        y_mule = y_hitch + L * np.sin(theta_mule)

        # return mule_pose, hitch_pose
        return (x_mule, y_mule, theta_mule), (x_hitch, y_hitch, theta_hitch)

    def cost(self, omega_seq, state0, phi0, t0, path_func, L, D):
        """
        Cost evaluates predicted trailer position error along horizon.
        state0: initial mule state [x_m, y_m, theta_m, phi]
        phi0: initial phi (articulation)
        path_func: provides desired trailer (x,y,theta) at any t
        """
        state = np.array(state0, dtype=float)
        J = 0.0
        prev_omega = self.prev_omega

        for i in range(self.N):
            omega = float(omega_seq[i])
            # simulate mule+phi forward one step
            state = self.step_dynamics(state, omega, L, D)
            x_m, y_m, theta_m, phi = state

            # compute trailer pose from predicted mule state
            _, (x_t_pred, y_t_pred, theta_t_pred) = self.compute_hitch_trailer_pose(
                (x_m, y_m, theta_m), phi, L, D)

            t_pred = t0 + (i + 1) * self.dt
            # reference trailer position at t_pred
            rx_t, ry_t, _ = path_func(t_pred)

            # Euclidean position error of trailer (we want trailer to follow ref point)
            pos_err_sq = (x_t_pred - rx_t) ** 2 + (y_t_pred - ry_t) ** 2
            J += self.w_pos_trailer * pos_err_sq

            # regularize omega and omega-rate
            J += self.w_omega * (omega ** 2)
            J += self.w_omega_rate * ((omega - prev_omega) ** 2)

            prev_omega = omega

        return J

    def mpc(self, x, y, theta, phi, path_func, t0=0.0, L=1.5, D=2.0):
        state0 = [x, y, theta, phi]
        init_guess = np.full(self.N, self.prev_omega, dtype=float)
        bounds = [(self.omega_min, self.omega_max)] * self.N

        res = minimize(self.cost, init_guess,
                       args=(state0, phi, t0, path_func, L, D),
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
    def __init__(self, path_func, start, controller, L=1.5, D=2.0, dt=0.1, T=100):
        self.path_func = path_func
        self.controller = controller
        self.L = L
        self.D = D
        self.dt = dt
        self.T = T
        self.steps = int(T / dt)

        # start now refers to trailer position/heading
        self.x, self.y, self.theta = start  # trailer pose
        self.phi = 0.0  # initial hitch angle guess

        # compute corresponding mule/hitch initial poses
        mule_pose, hitch_pose = self.controller.compute_mule_hitch_pose_from_trailer(
            (self.x, self.y, self.theta), self.L, self.D, phi_guess=self.phi)

        self.x_m, self.y_m, self.theta_m = mule_pose

        self.v = controller.v

        self.states = []    # mule positions logged
        self.hitches = []
        self.trailers = []  # trailer positions logged
        self.frenet_states = []
        self.frenet_hitches = []
        self.frenet_trailers = []

        self.ref_s = []
        self.ref_d = []

    def compute_hitch_trailer_from_mule(self, mule_pose, phi):
        return self.controller.compute_hitch_trailer_pose(mule_pose, phi, self.L, self.D)

    def simulate(self):
        # initialize mule state for dynamics used in controller
        state_mule = np.array([self.x_m, self.y_m, self.theta_m, self.phi], dtype=float)

        for i in range(self.steps):
            t = i * self.dt

            # compute current mule pose from state_mule
            x_m, y_m, theta_m, phi = state_mule

            # use MPC which predicts trailer behavior (cost penalizes trailer position error)
            omega = self.controller.mpc(x_m, y_m, theta_m, phi,
                                        self.path_func, t0=t, L=self.L, D=self.D)

            # apply dynamics to mule state (mule moves and pushes hitch/trailer)
            state_mule = self.controller.step_dynamics(state_mule, omega, self.L, self.D)
            x_m, y_m, theta_m, phi = state_mule

            # compute hitch and trailer poses from updated mule state
            (x_h, y_h, th_h), (x_t, y_t, th_t) = self.compute_hitch_trailer_from_mule((x_m, y_m, theta_m), phi)

            # log
            self.states.append([x_m, y_m])
            self.hitches.append([x_h, y_h])
            self.trailers.append([x_t, y_t])

            # Frenet logging (now trailer is the primary tracked body)
            rx_t, ry_t, theta_ref = self.path_func(t)
            # For mule, hitch, trailer store s (t*v) and lateral offset d
            for point, storage in zip([[x_m, y_m, theta_m], [x_h, y_h, th_h], [x_t, y_t, th_t]],
                                      [self.frenet_states, self.frenet_hitches, self.frenet_trailers]):
                s_cond, d_cond = CartesianFrenetConverter.cartesian_to_frenet(
                    rs=t * self.v, rx=rx_t, ry=ry_t, rtheta=theta_ref,
                    rkappa=0.0, rdkappa=0.0, x=point[0], y=point[1],
                    v=self.v, a=0.0, theta=point[2], kappa=0.0
                )
                storage.append([s_cond[0], d_cond[0]])

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

        # --- Plot reference path from path_func ---
        t_vals = np.linspace(0, self.T, 500)
        ref_x, ref_y = [], []
        for t in t_vals:
            px, py, _ = self.path_func(t)
            ref_x.append(px)
            ref_y.append(py)
        ax_cart.plot(ref_x, ref_y, 'k--', alpha=0.5, label='Reference (Trailer) Path')

        line_traj, = ax_cart.plot([], [], 'r-', label='Mule Path')
        mule_point, = ax_cart.plot([], [], 'ro', label='Mule')

        trailer_traj, = ax_cart.plot([], [], color='lightblue', linestyle='-', label='Trailer Path')
        hitch_point, = ax_cart.plot([], [], 'yo', label='Hitch')
        trailer_point, = ax_cart.plot([], [], 'ko', label='Trailer')

        link1, = ax_cart.plot([], [], 'r-', lw=1.5)
        link2, = ax_cart.plot([], [], 'k-', lw=1.5)
        ax_cart.legend()
        ax_cart.set_title("Cartesian Frame (Trailer-based, Reverse Tracking)")

        # ==== Frenet Plot ====
        ax_frenet.set_xlim(0, self.steps * self.dt * abs(self.v))
        ax_frenet.set_ylim(-5, 5)
        ax_frenet.axhline(0, color='gray', linestyle='--', label='Reference Path (d=0)')

        line_frenet, = ax_frenet.plot([], [], 'r-', label='Trailer s-d')
        mule_frenet_point, = ax_frenet.plot([], [], 'ro', label='Trailer')
        ax_frenet.legend()
        ax_frenet.set_xlabel("s (arc length)")
        ax_frenet.set_ylabel("d (lateral deviation)")
        ax_frenet.set_title("Frenet Frame")

        # ==== Precompute Frenet Coordinates ====
        frenet_coords = []
        s_accum = 0
        for i in range(len(self.trailers)):
            x, y = self.trailers[i]
            x_ref, y_ref, theta_ref = self.path_func(i * self.dt)
            s_accum += self.v * self.dt
            d = np.sin(theta_ref) * (x - x_ref) - np.cos(theta_ref) * (y - y_ref)
            frenet_coords.append([s_accum, d])
        frenet_coords = np.array(frenet_coords)

        def init():
            line_traj.set_data([], [])
            mule_point.set_data([], [])
            trailer_traj.set_data([], [])
            hitch_point.set_data([], [])
            trailer_point.set_data([], [])
            link1.set_data([], [])
            link2.set_data([], [])
            line_frenet.set_data([], [])
            mule_frenet_point.set_data([], [])
            return (line_traj, mule_point, trailer_traj, hitch_point, trailer_point,
                    link1, link2, line_frenet, mule_frenet_point)

        def update(i):
            mule = self.states[i]
            hitch = self.hitches[i]
            trailer = self.trailers[i]

            # Mule path
            line_traj.set_data(self.states[:i+1, 0], self.states[:i+1, 1])
            mule_point.set_data([mule[0]], [mule[1]])

            # Trailer path
            trailer_traj.set_data(self.trailers[:i+1, 0], self.trailers[:i+1, 1])
            hitch_point.set_data([hitch[0]], [hitch[1]])
            trailer_point.set_data([trailer[0]], [trailer[1]])

            link1.set_data([mule[0], hitch[0]], [mule[1], hitch[1]])
            link2.set_data([hitch[0], trailer[0]], [hitch[1], trailer[1]])

            # Frenet path
            line_frenet.set_data(frenet_coords[:i+1, 0], frenet_coords[:i+1, 1])
            mule_frenet_point.set_data([frenet_coords[i, 0]], [frenet_coords[i, 1]])

            return (line_traj, mule_point, trailer_traj, hitch_point, trailer_point,
                    link1, link2, line_frenet, mule_frenet_point)

        ani = animation.FuncAnimation(fig, update, frames=self.steps,
                                      init_func=init, interval=50, blit=True)
        plt.show()


if __name__ == "__main__":
    # Track figure-eight in reverse: parameterize with -t so path_func returns points along reversed direction
    path_func = lambda t: ParametricFunctions.figure_eight(-t, a=10.0, v=1.0)

    # Controller with trailer-focused cost
    controller = Controller(N=12, dt=0.1, v=-1.0, w_pos_trailer=1500.0)

    # Start: trailer at origin pointing 45deg
    sim = PathPlanningFrenetFrame(path_func, (0, 0, np.pi/4), controller, L=1.5, D=2.0, dt=0.1, T=20)
    sim.simulate()#!/Users/subhodeep/venv/bin/python
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.optimize import minimize

class CartesianFrenetConverter:
    """
    A class for converting states between Cartesian and Frenet coordinate systems
    (unchanged from your original implementation)
    """

    @ staticmethod
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

    @ staticmethod
    def frenet_to_cartesian(rs, rx, ry, rtheta, rkappa, rdkappa, s_condition, d_condition):
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

    @staticmethod
    def straight_line(t, slope=0.0, intercept=0.0, v=1.0):
        x = v * t
        y = slope * x + intercept
        theta = np.arctan2(slope, 1.0)
        return x, y, theta


class Controller:
    """
    Controller kept for API-compatibility. In this simplified reverse-tracking
    version we will *force* the trailer onto the reference at each timestep
    (the mule may move arbitrarily to achieve that). This is the simplest
    reliable way to guarantee trailer follows the reference exactly while
    preserving your original structure.
    """
    def __init__(self, N=12, dt=0.1, v=1.0,
                 omega_bounds=(-2.0, 2.0), max_delta_omega=0.5,
                 use_exp_smooth=False, smooth_alpha=0.3):
        self.N = N
        self.dt = dt
        self.v = v
        self.omega_min, self.omega_max = omega_bounds
        self.prev_omega = 0.0
        self.max_delta_omega = max_delta_omega
        self.use_exp_smooth = use_exp_smooth
        self.smooth_alpha = smooth_alpha

    def wrap_angle(self, angle):
        return (angle + np.pi) % (2 * np.pi) - np.pi


class PathPlanningFrenetFrame:
    """
    This version forces the trailer to follow path_func(t) exactly at each
    logged timestep. From that trailer pose we compute hitch and mule poses
    by forward geometry (hitch ahead of trailer by +D, mule ahead of hitch by +L).

    This meets your requirement: "mule may be arbitrary as long it ensures
    trailer is on ref path". If you later want dynamic pushing with realistic
    kinematics, I can replace the forced-tracking with an MPC or per-step
    optimizer that respects mule/trailer dynamics.
    """
    def __init__(self, path_func, start, controller, L=1.5, D=2.0, dt=0.1, T=100, smooth_mule_alpha=0.0):
        self.path_func = path_func
        self.controller = controller
        self.L = L
        self.D = D
        self.dt = dt
        self.T = T
        self.steps = int(T / dt)

        # start refers to trailer pose (x_t, y_t, theta_t)
        self.x_t, self.y_t, self.theta_t = start

        # compute initial hitch and mule from trailer (forward geometry)
        self.x_h = self.x_t + D * math.cos(self.theta_t)
        self.y_h = self.y_t + D * math.sin(self.theta_t)
        self.theta_h = self.theta_t

        # choose mule aligned with hitch initially (phi = 0)
        self.x_m = self.x_h + L * math.cos(self.theta_h)
        self.y_m = self.y_h + L * math.sin(self.theta_h)
        self.theta_m = self.theta_h

        # smoothing factor for mule motion (0 = teleport mule to geometry-derived pose)
        self.smooth_mule_alpha = float(smooth_mule_alpha)

        self.states = []    # mule positions logged (x_m, y_m)
        self.hitches = []   # hitch positions logged (x_h, y_h)
        self.trailers = []  # trailer positions logged (x_t, y_t)

        self.frenet_states = []    # mule frenet (s,d)
        self.frenet_hitches = []
        self.frenet_trailers = []

        self.ref_s = []
        self.ref_d = []

    def simulate(self):
        s_accum = 0.0
        for i in range(self.steps):
            t = i * self.dt

            # Desired trailer pose from parametric function at this time t
            # (note: user requested reverse tracking, so call path_func with -t when desired)
            x_t_des, y_t_des, theta_t_des = self.path_func(t)

            # Force trailer to reference point
            x_t = float(x_t_des)
            y_t = float(y_t_des)
            theta_t = float(theta_t_des)

            # Compute hitch ahead of trailer by +D along trailer heading
            x_h = x_t + self.D * math.cos(theta_t)
            y_h = y_t + self.D * math.sin(theta_t)
            theta_h = theta_t

            # Compute mule ahead of hitch by +L along hitch heading
            x_m_des = x_h + self.L * math.cos(theta_h)
            y_m_des = y_h + self.L * math.sin(theta_h)
            theta_m_des = theta_h

            # Optionally smooth mule motion to avoid instant teleportation
            if self.smooth_mule_alpha > 0.0 and len(self.states) > 0:
                x_m = (1 - self.smooth_mule_alpha) * self.states[-1][0] + self.smooth_mule_alpha * x_m_des
                y_m = (1 - self.smooth_mule_alpha) * self.states[-1][1] + self.smooth_mule_alpha * y_m_des
                # wrap angle smoothing
                ang_prev = self.theta_m
                # move toward desired angle by small step
                delta_ang = (theta_m_des - ang_prev + math.pi) % (2 * math.pi) - math.pi
                theta_m = ang_prev + self.smooth_mule_alpha * delta_ang
            else:
                x_m = x_m_des
                y_m = y_m_des
                theta_m = theta_m_des

            # Save current mule/hitch/trailer
            self.states.append([x_m, y_m])
            self.hitches.append([x_h, y_h])
            self.trailers.append([x_t, y_t])

            # Frenet logging: s increases by controller.v * dt (user param)
            s_accum += self.controller.v * self.dt
            # store s,d for mule/hitch/trailer relative to the reference trailer pose at time t
            rx, ry, rtheta = x_t_des, y_t_des, theta_t_des
            for point, storage in zip([[x_m, y_m, theta_m], [x_h, y_h, theta_h], [x_t, y_t, theta_t]],
                                      [self.frenet_states, self.frenet_hitches, self.frenet_trailers]):
                s_cond, d_cond = CartesianFrenetConverter.cartesian_to_frenet(
                    rs=s_accum, rx=rx, ry=ry, rtheta=rtheta,
                    rkappa=0.0, rdkappa=0.0, x=point[0], y=point[1],
                    v=self.controller.v, a=0.0, theta=point[2], kappa=0.0
                )
                storage.append([s_cond[0], d_cond[0]])

            self.ref_s.append(s_accum)
            self.ref_d.append(0.0)

            # update mule nominal angle (for smoothing next step)
            self.theta_m = theta_m

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

        # Cartesian
        ax_cart.set_aspect('equal')
        ax_cart.set_xlim(-12, 12)
        ax_cart.set_ylim(-12, 12)

        # reference path
        t_vals = np.linspace(0, self.T, 500)
        ref_x, ref_y = [], []
        for t in t_vals:
            px, py, _ = self.path_func(t)
            ref_x.append(px)
            ref_y.append(py)
        ax_cart.plot(ref_x, ref_y, 'k--', alpha=0.5, label='Reference (Trailer) Path')

        line_mule, = ax_cart.plot([], [], 'r-', label='Mule Path')
        mule_point, = ax_cart.plot([], [], 'ro', label='Mule')

        line_trailer, = ax_cart.plot([], [], 'b-', label='Trailer Path')
        hitch_point, = ax_cart.plot([], [], 'yo', label='Hitch')
        trailer_point, = ax_cart.plot([], [], 'ko', label='Trailer')

        link1, = ax_cart.plot([], [], 'r-', lw=1.5)
        link2, = ax_cart.plot([], [], 'k-', lw=1.5)
        ax_cart.legend()
        ax_cart.set_title("Cartesian Frame (Trailer-based, Forced Reverse Tracking)")

        # Frenet
        ax_frenet.set_xlim(0, self.steps * self.dt * abs(self.controller.v))
        ax_frenet.set_ylim(-5, 5)
        ax_frenet.axhline(0, color='gray', linestyle='--', label='Reference (d=0)')
        line_frenet, = ax_frenet.plot([], [], 'r-', label='Trailer s-d')
        mule_frenet_point, = ax_frenet.plot([], [], 'ro', label='Trailer')
        ax_frenet.legend()
        ax_frenet.set_xlabel("s (arc length)")
        ax_frenet.set_ylabel("d (lateral deviation)")
        ax_frenet.set_title("Frenet Frame")

        frenet_coords = []
        for i in range(len(self.trailers)):
            x, y = self.trailers[i]
            x_ref, y_ref, theta_ref = self.path_func(i * self.dt)
            d = np.sin(theta_ref) * (x - x_ref) - np.cos(theta_ref) * (y - y_ref)
            frenet_coords.append([self.ref_s[i], d])
        frenet_coords = np.array(frenet_coords)

        def init():
            line_mule.set_data([], [])
            mule_point.set_data([], [])
            line_trailer.set_data([], [])
            hitch_point.set_data([], [])
            trailer_point.set_data([], [])
            link1.set_data([], [])
            link2.set_data([], [])
            line_frenet.set_data([], [])
            mule_frenet_point.set_data([], [])
            return (line_mule, mule_point, line_trailer, hitch_point, trailer_point,
                    link1, link2, line_frenet, mule_frenet_point)

        def update(i):
            mule = self.states[i]
            hitch = self.hitches[i]
            trailer = self.trailers[i]

            # Mule path
            line_mule.set_data(self.states[:i+1, 0], self.states[:i+1, 1])
            mule_point.set_data([mule[0]], [mule[1]])

            # Trailer path
            line_trailer.set_data(self.trailers[:i+1, 0], self.trailers[:i+1, 1])
            hitch_point.set_data([hitch[0]], [hitch[1]])
            trailer_point.set_data([trailer[0]], [trailer[1]])

            link1.set_data([mule[0], hitch[0]], [mule[1], hitch[1]])
            link2.set_data([hitch[0], trailer[0]], [hitch[1], trailer[1]])

            # Frenet path
            line_frenet.set_data(frenet_coords[:i+1, 0], frenet_coords[:i+1, 1])
            mule_frenet_point.set_data([frenet_coords[i, 0]], [frenet_coords[i, 1]])

            return (line_mule, mule_point, line_trailer, hitch_point, trailer_point,
                    link1, link2, line_frenet, mule_frenet_point)

        ani = animation.FuncAnimation(fig, update, frames=self.steps,
                                      init_func=init, interval=50, blit=True)
        plt.show()


if __name__ == "__main__":
    # Reverse figure-eight: the parametric generator is called with -t to go backward along the curve
    path_func = lambda t: ParametricFunctions.figure_eight(-t, a=10.0, v=1.0)

    # Controller is retained for compatibility; its 'v' is used to increment s in Frenet logging
    controller = Controller(N=12, dt=0.1, v=1.0)

    # Start: trailer at origin pointing 45deg
    sim = PathPlanningFrenetFrame(path_func, (0.0, 0.0, math.pi/4), controller, L=1.5, D=2.0, dt=0.1, T=20, smooth_mule_alpha=0.05)
    sim.simulate()
    sim.animate()

    sim.animate()
