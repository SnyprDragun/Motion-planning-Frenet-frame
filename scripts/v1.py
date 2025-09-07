#!/Users/subhodeep/venv/bin/python

import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from dataclasses import dataclass
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

        cos_theta_r = np.cos(rtheta)
        sin_theta_r = np.sin(rtheta)

        cross_rd_nd = cos_theta_r * dy - sin_theta_r * dx
        d = np.copysign(np.hypot(dx, dy), cross_rd_nd)

        delta_theta = theta - rtheta
        tan_delta_theta = np.tan(delta_theta)
        cos_delta_theta = np.cos(delta_theta)

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

        cos_theta_r = np.cos(rtheta)
        sin_theta_r = np.sin(rtheta)

        x = rx - sin_theta_r * d_condition[0]
        y = ry + cos_theta_r * d_condition[0]

        one_minus_kappa_r_d = 1 - rkappa * d_condition[0]

        tan_delta_theta = d_condition[1] / one_minus_kappa_r_d
        delta_theta = np.atan2(d_condition[1], one_minus_kappa_r_d)
        cos_delta_theta = np.cos(delta_theta)

        theta = CartesianFrenetConverter.normalize_angle(delta_theta + rtheta)

        kappa_r_d_prime = rdkappa * d_condition[0] + rkappa * d_condition[1]

        kappa = (((d_condition[2] + kappa_r_d_prime * tan_delta_theta) *
                  cos_delta_theta * cos_delta_theta) / one_minus_kappa_r_d + rkappa) * \
            cos_delta_theta / one_minus_kappa_r_d

        d_dot = d_condition[1] * s_condition[1]
        v = np.sqrt(one_minus_kappa_r_d * one_minus_kappa_r_d *
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
        return (angle + np.pi) % (2 * np.pi) - np.pi


class ParametricFunctions:
    """
    Generates parametric paths. Now also returns the path's curvature.
    """
    @staticmethod
    def circle(t, R=10.0, v=1.0):
        """
        Generates a circle path and its kinematic properties.
        Returns: x, y, theta (heading), kappa (curvature)
        """
        omega = v / R
        x = R * np.cos(omega * t)
        y = R * np.sin(omega * t)
        theta = np.arctan2(y, x) + np.pi / 2
        kappa = 1.0 / R
        return x, y, theta, kappa

    @staticmethod
    def ellipse(t, a=12.0, b=8.0, v=1.0):
        """
        Generates an ellipse path and its kinematic properties.
        Returns: x, y, theta (heading), kappa (curvature)
        """
        # Approximate constant velocity
        omega = v / ((a + b) / 2.0)
        
        # Position
        x = a * np.cos(omega * t)
        y = b * np.sin(omega * t)
        
        # First derivatives
        dx = -a * omega * np.sin(omega * t)
        dy = b * omega * np.cos(omega * t)
        
        # Second derivatives
        ddx = -a * omega**2 * np.cos(omega * t)
        ddy = -b * omega**2 * np.sin(omega * t)
        
        # Heading
        theta = np.arctan2(dy, dx)
        
        # Curvature
        speed_sq = dx**2 + dy**2
        if speed_sq < 1e-6:
            kappa = 0
        else:
            kappa = (dx * ddy - dy * ddx) / speed_sq**(1.5)
            
        return x, y, theta, kappa

    @staticmethod
    def figure_eight(t, a=10.0, v=1.0):
        """
        Generates a figure-eight path and its kinematic properties.
        Returns: x, y, theta (heading), kappa (curvature)
        """
        omega = v / a
        
        # First derivatives (velocity components)
        dx = a * omega * np.cos(omega * t)
        dy = a * omega * np.cos(2 * omega * t)
        
        # Second derivatives (acceleration components)
        ddx = -a * omega**2 * np.sin(omega * t)
        ddy = -2 * a * omega**2 * np.sin(2 * omega * t)

        # Path properties
        x = a * np.sin(omega * t)
        y = a * np.sin(omega * t) * np.cos(omega * t)
        theta = np.arctan2(dy, dx)
        
        speed_sq = dx**2 + dy**2
        if speed_sq < 1e-6:
            kappa = 0
        else:
            kappa = (dx * ddy - dy * ddx) / speed_sq**(1.5)

        return x, y, theta, kappa

    @staticmethod
    def straight_line(t, slope=0.0, intercept=0.0, v=1.0):
        """
        Generates a straight line path and its kinematic properties.
        Returns: x, y, theta (heading), kappa (curvature)
        """
        x = v * t
        y = slope * x + intercept
        theta = np.arctan2(slope, 1.0)
        kappa = 0.0
        return x, y, theta, kappa


class Controller:
    """
    MPC Controller with fine-tuned weights for smooth and stable motion.
    """
    def __init__(self, N=12, dt=0.1, v=1.0,
                 w_d_mule=1200.0, w_d_hitch=0.0, w_d_trailer=0.0, 
                 w_d_dot_trailer=200.0, w_phi=800.0,
                 w_omega_mule=0.01, w_omega_rate_mule=1.0,
                 omega_bounds=(-2.0, 2.0), max_delta_omega=0.5,
                 use_exp_smooth=False, smooth_alpha=0.3, is_reverse=False):
        """
        MPC controller using only lateral Frenet errors of mule, hitch, trailer.
        - N: horizon steps
        - dt: timestep
        - lookahead = N * dt
        - v: forward speed (constant)
        - w_d_mule/hitch/trailer: weights on lateral deviations
        - w_d_mule is the key weight to enforce strict lateral following (increase if you want stricter)
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
        self.w_d_dot_trailer = w_d_dot_trailer
        self.w_phi = w_phi

        self.w_omega_mule = w_omega_mule
        self.w_omega_rate_mule = w_omega_rate_mule

        self.omega_min, self.omega_max = omega_bounds
        self.prev_omega = 0.0
        self.max_delta_omega = max_delta_omega

        self.use_exp_smooth = use_exp_smooth
        self.smooth_alpha = smooth_alpha

        self.is_reverse = is_reverse

        if self.is_reverse:
            self.w_d_trailer = 3500.0
            self.w_omega_rate_mule = 9000.0
        else:
            self.w_d_trailer = 0.0
            self.w_omega_rate_mule=1.0

    def step_dynamics_forward(self, state, omega, L, D):
        """
        Mule frame
        Propagates the mule-centric state forward in time for one step.
        State: [x_mule, y_mule, theta_mule, phi_hitch]
        """
        x_mule, y_mule, theta_mule, phi_hitch = state
        dx = self.v * np.cos(theta_mule)
        dy = self.v * np.sin(theta_mule)
        dtheta = omega
        dphi = (self.v * np.sin(theta_mule - phi_hitch) - L * omega * np.cos(theta_mule - phi_hitch)) / D

        x_mule_new = x_mule + dx * self.dt
        y_mule_new = y_mule + dy * self.dt
        theta_mule_new = CartesianFrenetConverter.normalize_angle(theta_mule + dtheta * self.dt)
        phi_hitch_new = CartesianFrenetConverter.normalize_angle(phi_hitch + dphi * self.dt)

        return np.array([
            x_mule_new,
            y_mule_new,
            theta_mule_new,
            phi_hitch_new
        ])

    def step_dynamics_reverse(self, state, omega, L, D):
        """
        Trailer frame
        Propagates the trailer-centric state forward in time for one step.
        State: [x_trailer, y_trailer, theta_trailer, phi_hitch]
        """
        x_trailer, y_trailer, theta_trailer, phi_hitch = state
        theta_mule = CartesianFrenetConverter.normalize_angle(theta_trailer + phi_hitch)
        vx_h = -self.v * np.cos(theta_mule) + L * omega * np.sin(theta_mule)
        vy_h = -self.v * np.sin(theta_mule) - L * omega * np.cos(theta_mule)
        v_trailer = vx_h * np.cos(theta_trailer) + vy_h * np.sin(theta_trailer)
        v_perp = -vx_h * np.sin(theta_trailer) + vy_h * np.cos(theta_trailer)
        dtheta_trailer = v_perp / D
        dphi = omega - dtheta_trailer

        x_trailer_new = x_trailer + v_trailer * np.cos(theta_trailer) * self.dt
        y_trailer_new = y_trailer + v_trailer * np.sin(theta_trailer) * self.dt
        theta_trailer_new = CartesianFrenetConverter.normalize_angle(theta_trailer + dtheta_trailer * self.dt)
        phi_hitch_new = CartesianFrenetConverter.normalize_angle(phi_hitch + dphi * self.dt)

        return np.array([
            x_trailer_new, 
            y_trailer_new, 
            theta_trailer_new, 
            phi_hitch_new])

    def lateral_deviation(self, x, y, theta_ref, rx, ry):
        """
        Compute lateral deviation d for pose (x,y) relative to reference path point (rx, ry, theta_ref).
        """
        return np.sin(theta_ref) * (x - rx) - np.cos(theta_ref) * (y - ry)

    def compute_hitch_trailer_pose(self, mule_pose, phi, L, D):
        """
        Compute hitch and trailer poses from mule pose and hitch angle phi.
        Returns: (x_hitch, y_hitch, theta_hitch), (x_trailer, y_trailer, theta_trailer)
        """
        x, y, theta = mule_pose

        # Hitch position behind mule by L
        x_hitch = x - L * np.cos(theta)
        y_hitch = y - L * np.sin(theta)
        theta_hitch = CartesianFrenetConverter.normalize_angle(theta + phi)

        # Trailer position behind hitch by D
        x_trailer = x_hitch - D * np.cos(theta_hitch)
        y_trailer = y_hitch - D * np.sin(theta_hitch)
        theta_trailer = theta_hitch  # trailer heading same as hitch

        return (x_hitch, y_hitch, theta_hitch), (x_trailer, y_trailer, theta_trailer)

    def cost(self, omega_seq, state0, t0, path_func, L, D):
        """
        Cost function using Frenet frame errors and hard constraints for robust tracking.
        """
        state = np.array(state0, dtype=float)
        J = 0.0
        prev_omega = self.prev_omega

        for i in range(self.N):
            omega = float(omega_seq[i])
            x_lead, y_lead, theta_lead, phi_hitch = state
            t_pred = t0 + (i + 1) * self.dt
            x_ref, y_ref, theta_ref, kappa_ref = path_func(t_pred)

            if self.is_reverse:
                # =============== for reverse ===============
                state = self.step_dynamics_reverse(state, omega, L, D)
                x_ref_trailer, y_ref_trailer, theta_ref_trailer, phi_hitch = state

                s_cond, d_cond = CartesianFrenetConverter.cartesian_to_frenet(
                    rs=t_pred * self.v, rx=x_ref, ry=y_ref, rtheta=theta_ref,
                    rkappa=kappa_ref, rdkappa=0.0, x=x_ref_trailer, y=y_ref_trailer,
                    v=self.v, a=0.0, theta=theta_ref_trailer, kappa=0.0
                )
                
                if d_cond is None:
                    J += 1e7
                    continue

                # Cost on lateral deviation (d) and lateral velocity (d_prime or d_dot).
                # Penalizing d_prime implicitly corrects heading error in a stable way.
                J += self.w_d_trailer * (d_cond[0] ** 2)
                J += self.w_d_dot_trailer * (d_cond[1] ** 2)
                
                # Soft cost on hitch angle
                J += self.w_phi * (phi_hitch ** 2)
                
                # *** CRITICAL: Hard constraint to prevent jack-knifing ***
                # Apply an enormous penalty if the hitch angle exceeds a safe physical limit.
                if abs(phi_hitch) > np.radians(85): # 85 degrees is a safe limit
                    J += 1e8 

                # Regularize control inputs for smoothness
                J += self.w_omega_mule * (omega ** 2)
                J += self.w_omega_rate_mule * ((omega - prev_omega) ** 2)
                prev_omega = omega

            else:
                # =============== for forward ===============
                state = self.step_dynamics_forward(state, omega, L, D)

                # Time-shifted reference points for hitch and trailer (accounting for physical offset)
                t_pred_hitch = max(t_pred - L / self.v, 0.0)
                t_pred_trailer = max(t_pred - (L + D) / self.v, 0.0)

                x_ref_hitch, y_ref_hitch, theta_ref_hitch, kappa_ref_hitch = path_func(t_pred_hitch)
                x_ref_trailer, y_ref_trailer, theta_ref_trailer, kappa_ref_trailer = path_func(t_pred_trailer)
                (x_hitch, y_hitch, theta_hitch), (x_trailer, y_trailer, theta_trailer) = self.compute_hitch_trailer_pose((x_lead, y_lead, theta_lead), phi_hitch, L, D)

                d_mule = self.lateral_deviation(x_lead, y_lead, theta_ref, x_ref, y_ref)
                d_hitch = self.lateral_deviation(x_hitch, y_hitch, theta_ref_hitch, x_ref_hitch, y_ref_hitch)
                d_trailer = self.lateral_deviation(x_trailer, y_trailer, theta_ref_trailer, x_ref_trailer, y_ref_trailer)

                # Accumulate weighted lateral deviation costs
                J += self.w_d_mule * (d_mule ** 2)
                J += self.w_d_hitch * (d_hitch ** 2)
                J += self.w_d_trailer * (d_trailer ** 2)

                # Regularize omega and omega rate change
                J += self.w_omega_mule * (omega ** 2)
                J += self.w_omega_rate_mule * ((omega - prev_omega) ** 2)

                heading_err = np.atan2(np.sin(theta_lead - theta_ref), np.cos(theta_lead - theta_ref))
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


@dataclass
class Obstacle:
    x: float
    y: float
    radius: float


class PathPlanningFrenetFrame:
    """
    Quintic-framed detour:
        - uses a quintic smoothstep phi(z)=10z^3-15z^4+6z^5 mapped to d(s)=A*phi(z)
        - pre-starts detour by `pre_start` meters so vehicle begins steering earlier
        - chooses sign to move away from obstacle and clamps amplitude
        - ensures d, d', d'' == 0 at s0 and s1 for smooth curvature continuity
    """
    def __init__(self, path_func, start, controller, L=1.5, D=2.0, dt=0.1, T=100, is_reverse=False):
        self.start_simulation = time.time()
        self.path_func = path_func
        self.controller = controller
        self.L = L
        self.D = D
        self.dt = dt
        self.T = T
        self.steps = int(T / dt)
        self.v = controller.v
        self.start = start
        self.x, self.y, self.theta = self.start
        self.is_reverse = is_reverse

        # obstacles + avoidance params
        self.obstacles = []
        self.clearance = 1.0
        self.influence_margin = 3.0
        self.detour_span = 6.0      # meters of s coverage (increase for smoother gentler detour)
        self.pre_start = 2.0        # meters *before* obstacle center where detour should begin
        self.max_amplitude = 2.5    # clamp amplitude to avoid extreme maneuvers (m)

        # detour state (in s)
        self.detour_active = False
        self.detour_s0 = None
        self.detour_s1 = None
        self.detour_s_center = None
        self.detour_A = None
        self.detour_obs = None

        self.sim_time = 0.0
        self.debug = False

        # logs
        self.mules = []
        self.hitches = []
        self.trailers = []
        self.frenet_mules = []
        self.frenet_hitches = []
        self.frenet_trailers = []

        if self.is_reverse == True:
            self.phi = 0.0

            # Pre-compute a high-resolution reference path for accurate plotting
            self.ref_path_dense = None
            self.ref_s_dense = None
            self._precompute_ref_path()

        else:
            self.phi = self.theta

            # Precompute Frenet reference path
            self.ref_s = []
            self.ref_d = []

    def add_obstacle(self, obstacle: Obstacle, clearance: float = 1.0, influence_margin: float = 3.0,
                     detour_span: float = 6.0, pre_start: float = 2.0, max_amplitude: float = 2.5):
        self.obstacles.append(obstacle)
        self.clearance = clearance
        self.influence_margin = influence_margin
        self.detour_span = detour_span
        self.pre_start = pre_start
        self.max_amplitude = max_amplitude

    def _precompute_ref_path(self):
        """Creates a high-resolution lookup table of the reference path."""
        # Create a dense set of points along the path, e.g., 5 points per simulation step
        num_dense_points = self.steps * 5
        t_ref = np.linspace(0, self.T, num_dense_points)

        self.ref_path_dense = np.array([self.path_func(ti) for ti in t_ref])

        # Calculate the true cumulative arc length 's' along this dense path
        path_diffs = np.diff(self.ref_path_dense[:, :2], axis=0)
        segment_lengths = np.sqrt(np.sum(path_diffs**2, axis=1))
        self.ref_s_dense = np.insert(np.cumsum(segment_lengths), 0, 0)

    def compute_mule_hitch_from_trailer(self):
        """Helper to compute mule/hitch poses from the current trailer state."""
        theta_m = CartesianFrenetConverter.normalize_angle(self.theta + self.phi)
        x_h = self.x + self.D * np.cos(self.theta)
        y_h = self.y + self.D * np.sin(self.theta)
        x_m = x_h + self.L * np.cos(theta_m)
        y_m = y_h + self.L * np.sin(theta_m)
        return np.array([x_m, y_m]), np.array([x_h, y_h])

    def compute_hitch_trailer(self):
        x_h = self.x - self.L * np.cos(self.theta)
        y_h = self.y - self.L * np.sin(self.theta)
        x_t = x_h - self.D * np.cos(self.phi)
        y_t = y_h - self.D * np.sin(self.phi)
        return np.array([x_h, y_h]), np.array([x_t, y_t])

    def _s_of_time(self, t):
        return t * self.v

    def _path_point(self, t):
        return self.path_func(t)

    def _distance_point_to_obs(self, x, y, obs):
        return np.hypot(x - obs.x, y - obs.y)

    def _should_plan_detour(self, current_t):
        lookahead = self.controller.N * self.controller.dt
        steps = max(1, int(lookahead / self.dt))
        for k in range(steps):
            t_check = current_t + k * self.dt
            px, py, ptheta, pkappa = self._path_point(t_check)
            for obs in self.obstacles:
                inner = obs.radius + self.clearance
                outer = inner + self.influence_margin
                d = self._distance_point_to_obs(px, py, obs)
                if d <= outer:
                    if self.debug:
                        print(f"[detect] t={current_t:.2f} found path near obs at t_check={t_check:.2f} d={d:.2f}")
                    return obs, t_check
        return None, None

    def _find_s_center(self, detect_t, obs, search_span=6.0):
        t0 = max(0.0, detect_t - search_span)
        t1 = detect_t + search_span
        N = max(20, int((t1 - t0) / self.dt))
        best_t = detect_t
        best_d = float('inf')
        for i in range(N + 1):
            tt = t0 + (t1 - t0) * (i / N)
            px, py, ptheta, pkappa = self._path_point(tt)
            d = self._distance_point_to_obs(px, py, obs)
            if d < best_d:
                best_d = d
                best_t = tt
        s_center = self._s_of_time(best_t)
        if self.debug:
            print(f"[center] best_t={best_t:.2f} best_d={best_d:.2f} s_center={s_center:.2f}")
        return s_center, best_d

    def _plan_quintic_detour(self, detect_t, obs):
        """
        Plan a quintic smoothstep detour in s-domain:
          - compute s_center near obstacle
          - set s0 = s_center - span/2 - pre_start (start earlier)
          - set s1 = s_center + span/2
          - amplitude A = sign * (obs.radius + clearance + trailer_margin) clamped to max_amplitude
          - use phi(z)=10z^3 -15z^4 +6z^5 so phi'(0)=phi'(1)=phi''(0)=phi''(1)=0
        """
        s_center, best_d = self._find_s_center(detect_t, obs)
        half = 0.5 * self.detour_span
        s0 = max(0.0, s_center - half - self.pre_start)
        s1 = s_center + half

        # compute sign: if obstacle is left of path, detour to right (negative d), else left.
        t_center = s_center / self.v
        rx, ry, theta_ref, kappa_ref = self._path_point(t_center)
        vx = obs.x - rx
        vy = obs.y - ry
        nx = -np.sin(theta_ref)
        ny = np.cos(theta_ref)
        dot = vx * nx + vy * ny
        sign = -1.0 if dot > 0.0 else 1.0

        # extra margin to account for trailer/hitch: project trailer endpoint clearance approx
        trailer_margin = max(0.0, self.D * 0.1 + self.L * 0.05)  # small heuristic

        A_req = sign * (obs.radius + self.clearance + trailer_margin)
        # clamp amplitude
        if abs(A_req) > self.max_amplitude:
            A = np.copysign(self.max_amplitude, A_req)
        else:
            A = A_req

        # store
        self.detour_active = True
        self.detour_s0 = s0
        self.detour_s1 = s1
        self.detour_s_center = s_center
        self.detour_A = A
        self.detour_obs = obs
        if self.debug:
            print(f"[plan_quintic] s0={s0:.2f} s_center={s_center:.2f} s1={s1:.2f} A={A:.2f} sign={sign}")

    def _maybe_deactivate_detour(self, current_s):
        if not self.detour_active:
            return
        if current_s > self.detour_s1 + 1e-2:
            if self.debug:
                print(f"[done] current_s={current_s:.2f} > s1={self.detour_s1:.2f} -> deactivate")
            self.detour_active = False
            self.detour_s0 = None
            self.detour_s1 = None
            self.detour_s_center = None
            self.detour_A = None
            self.detour_obs = None

    # quintic smoothstep phi and derivative helpers
    @staticmethod
    def _phi(z):
        # phi(z) = 10 z^3 - 15 z^4 + 6 z^5
        return 10*z**3 - 15*z**4 + 6*z**5

    @staticmethod
    def _phi_d(z):
        # phi'(z) = 30 z^2 - 60 z^3 + 30 z^4
        return 30*z**2 - 60*z**3 + 30*z**4

    @staticmethod
    def _phi_dd(z):
        # phi''(z) = 60 z - 180 z^2 + 120 z^3
        return 60*z - 180*z**2 + 120*z**3

    def _quintic_bump(self, s):
        """
        d(s) = A * phi(z), z=(s - s0)/(s1 - s0)
        returns d, d_s, d_ss
        d_s = A * phi'(z) * (1 / L)
        d_ss = A * phi''(z) * (1 / L^2)
        """
        if (not self.detour_active) or (self.detour_s0 is None):
            return 0.0, 0.0, 0.0
        s0 = self.detour_s0
        s1 = self.detour_s1
        if s < s0 - 1e-9 or s > s1 + 1e-9:
            return 0.0, 0.0, 0.0
        L = max(1e-6, s1 - s0)
        z = (s - s0) / L
        A = self.detour_A
        phi = self._phi(z)
        phi_d = self._phi_d(z)
        phi_dd = self._phi_dd(z)
        d = A * phi
        d_s = A * phi_d / L
        d_ss = A * phi_dd / (L * L)
        return d, d_s, d_ss

    def avoidance_path_frenet(self, t_pred):
        rx, ry, rtheta, rkappa = self._path_point(t_pred)
        s = self._s_of_time(t_pred)
        if not self.detour_active:
            return rx, ry, rtheta, rkappa
        if s < (self.detour_s0 - 1e-9) or s > (self.detour_s1 + 1e-9):
            return rx, ry, rtheta, rkappa
        d, d_s, d_ss = self._quintic_bump(s)
        s_cond = [s, self.v, 0.0]
        d_cond = [d, d_s, d_ss]
        try:
            x, y, theta_det, kappa, v_out, a_out = CartesianFrenetConverter.frenet_to_cartesian(
                rs=s, rx=rx, ry=ry, rtheta=rtheta, rkappa=0.0, rdkappa=0.0,
                s_condition=s_cond, d_condition=d_cond)
        except Exception as e:
            if self.debug:
                print(f"[conv_err] {e}")
            return rx, ry, rtheta, rkappa
        return x, y, theta_det, kappa

    def simulate(self):
        if self.is_reverse == True:
            last_closest_idx = 0
            for i in range(self.steps):
                t = i * self.dt
                
                omega = self.controller.mpc(self.x, self.y, self.theta, self.phi,
                                            self.path_func, t0=t, L=self.L, D=self.D)

                # Update state using reverse dynamics
                theta_m = CartesianFrenetConverter.normalize_angle(self.theta + self.phi)
                vx_h = -self.v * np.cos(theta_m) + self.L * omega * np.sin(theta_m)
                vy_h = -self.v * np.sin(theta_m) - self.L * omega * np.cos(theta_m)
                v_t = vx_h * np.cos(self.theta) + vy_h * np.sin(self.theta)
                v_perp = -vx_h * np.sin(self.theta) + vy_h * np.cos(self.theta)
                dtheta_dt = v_perp / self.D
                dphi_dt = omega - dtheta_dt
                
                self.x += v_t * np.cos(self.theta) * self.dt
                self.y += v_t * np.sin(self.theta) * self.dt
                self.theta = CartesianFrenetConverter.normalize_angle(self.theta + dtheta_dt * self.dt)
                self.phi = CartesianFrenetConverter.normalize_angle(self.phi + dphi_dt * self.dt)

                # Log Cartesian poses for plotting
                mule, hitch = self.compute_mule_hitch_from_trailer()
                self.trailers.append([self.x, self.y])
                self.hitches.append(hitch)
                self.mules.append(mule)

                # *** ROBUST FRENET PLOTTING LOGIC ***
                # Find the geometrically closest point on the pre-computed reference path
                current_pos = np.array([self.x, self.y])
                
                # Search in a window around the last known closest point for efficiency
                search_radius = 100
                start_idx = max(0, last_closest_idx - search_radius)
                end_idx = min(len(self.ref_path_dense), last_closest_idx + search_radius)
                
                distances_sq = np.sum((self.ref_path_dense[start_idx:end_idx, :2] - current_pos)**2, axis=1)
                
                # Find the index of the minimum distance in the search window
                local_closest_idx = np.argmin(distances_sq)
                closest_idx = start_idx + local_closest_idx
                last_closest_idx = closest_idx

                # Get the properties of the closest reference point
                closest_ref_point = self.ref_path_dense[closest_idx]
                rx, ry, rtheta = closest_ref_point[0], closest_ref_point[1], closest_ref_point[2]
                s_true = self.ref_s_dense[closest_idx]

                # Calculate lateral deviation 'd' relative to the closest point
                dx = self.x - rx
                dy = self.y - ry
                d = np.copysign(np.hypot(dx, dy), np.sin(rtheta) * dx - np.cos(rtheta) * dy)
                
                self.frenet_trailers.append([s_true, d])


            # Convert lists to numpy arrays
            self.trailers = np.array(self.trailers)
            self.hitches = np.array(self.hitches)
            self.mules = np.array(self.mules)
            self.frenet_trailers = np.array(self.frenet_trailers)

        else:
            for i in range(self.steps):
                t = i * self.dt

                current_s = self._s_of_time(t)

                if not self.detour_active:
                    obs, detect_t = self._should_plan_detour(t)
                    if obs is not None and detect_t is not None:
                        self._plan_quintic_detour(detect_t, obs)

                self._maybe_deactivate_detour(current_s)

                path_func_local = lambda t_pred: self.avoidance_path_frenet(t_pred)

                omega = self.controller.mpc(self.x, self.y, self.theta, self.phi,
                                            path_func_local, t0=t, L=self.L, D=self.D)

                # dynamics
                self.phi += ((self.v * np.sin(self.theta - self.phi) - self.L * omega * np.cos(self.theta - self.phi)) / self.D) * self.dt
                self.x += self.v * np.cos(self.theta) * self.dt
                self.y += self.v * np.sin(self.theta) * self.dt
                self.theta += omega * self.dt
                self.theta = CartesianFrenetConverter.normalize_angle(self.theta)

                self.sim_time += self.dt

                hitch, trailer = self.compute_hitch_trailer()

                # Log Cartesian
                self.mules.append([self.x, self.y])
                self.hitches.append(hitch)
                self.trailers.append(trailer)

                # Log Frenet (for mule, hitch, trailer) using original parametric ref
                x_ref, y_ref, theta_ref, kappa_ref = self.path_func(t)
                for point, storage in zip([[self.x, self.y, self.theta],
                                        [hitch[0], hitch[1], self.theta],
                                        [trailer[0], trailer[1], self.theta]],
                                        [self.frenet_mules, self.frenet_hitches, self.frenet_trailers]):
                    s_cond, d_cond = CartesianFrenetConverter.cartesian_to_frenet(
                        rs=t*self.v, rx=x_ref, ry=y_ref, rtheta=theta_ref,
                        rkappa=0.0, rdkappa=0.0, x=point[0], y=point[1],
                        v=self.v, a=0.0, theta=point[2], kappa=0.0
                    )
                    storage.append([s_cond[0], d_cond[0]])

                self.ref_s.append(current_s)
                self.ref_d.append(0.0)

            # convert to arrays
            self.mules = np.array(self.mules)
            self.hitches = np.array(self.hitches)
            self.trailers = np.array(self.trailers)
            self.frenet_mules = np.array(self.frenet_mules)
            self.frenet_hitches = np.array(self.frenet_hitches)
            self.frenet_trailers = np.array(self.frenet_trailers)
            self.ref_s = np.array(self.ref_s)
            self.ref_d = np.array(self.ref_d)

    def animate(self):
        fig = plt.figure(figsize=(16, 8))
        ax_cartesian = fig.add_subplot(1, 2, 1)
        ax_frenet = fig.add_subplot(1, 2, 2)

        # ==== Cartesian Plot ====
        ax_cartesian.set_aspect('equal')
        ax_cartesian.set_xlim(-12, 12)
        ax_cartesian.set_ylim(-12, 12)

        # --- Plot reference path from path_func ---
        t_vals = np.linspace(0, self.T, 500)
        ref_x, ref_y = [], []
        for t in t_vals:
            px, py, ptheta, pkappa = self.path_func(t)
            ref_x.append(px)
            ref_y.append(py)
        ax_cartesian.plot(ref_x, ref_y, 'k--', alpha=0.5, label='Reference Path')

        if self.is_reverse == True:
            line_trailer, = ax_cartesian.plot([], [], color='b', linestyle='-', label='Trailer Path')
            trailer_point, = ax_cartesian.plot([], [], 'ko', ms=8, label='Trailer (Lead)')
            hitch_point, = ax_cartesian.plot([], [], 'yo', ms=6, label='Hitch')
            mule_point, = ax_cartesian.plot([], [], 'ro', ms=8, label='Mule (Steering)')
            link1, = ax_cartesian.plot([], [], 'k-', lw=2)
            link2, = ax_cartesian.plot([], [], 'k-', lw=2)
            ax_cartesian.legend(loc='upper right')
            ax_cartesian.set_title("Reverse Path Tracking (Cartesian Frame)")
            ax_cartesian.grid(True)

            # === Frenet Error Plot ===
            max_s = self.ref_s_dense[-1]
            ax_frenet.set_xlim(0, max_s)
            ax_frenet.set_ylim(-5, 5) # Reasonable error bounds
            ax_frenet.axhline(0, color='gray', linestyle='--')
            frenet_trailer_error, = ax_frenet.plot([], [], 'g-', label='Lateral Error (d)')
            ax_frenet.set_title("Frenet Frame Lateral Error")
            ax_frenet.set_xlabel("True Arc Length (s) [m]")
            ax_frenet.set_ylabel("Lateral Deviation (d) [m]")
            ax_frenet.legend()
            ax_frenet.grid(True)

            fig.tight_layout()

            def init():
                line_trailer.set_data([], [])
                trailer_point.set_data([], [])
                hitch_point.set_data([], [])
                mule_point.set_data([], [])
                link1.set_data([], [])
                link2.set_data([], [])
                frenet_trailer_error.set_data([], [])
                return (line_trailer, trailer_point, hitch_point, mule_point, link1, link2, frenet_trailer_error)

            def update(i):
                trailer = self.trailers[i]
                hitch = self.hitches[i]
                mule = self.mules[i]

                line_trailer.set_data(self.trailers[:i+1, 0], self.trailers[:i+1, 1])
                trailer_point.set_data([trailer[0]], [trailer[1]])
                hitch_point.set_data([hitch[0]], [hitch[1]])
                mule_point.set_data([mule[0]], [mule[1]])
                link1.set_data([trailer[0], hitch[0]], [trailer[1], hitch[1]])
                link2.set_data([hitch[0], mule[0]], [hitch[1], mule[1]])

                # Corrected Frenet plot update
                if len(self.frenet_trailers) > i:
                    frenet_trailer_data = self.frenet_trailers[:i+1]
                    frenet_trailer_error.set_data(frenet_trailer_data[:, 0], frenet_trailer_data[:, 1])

                return (line_trailer, trailer_point, hitch_point, mule_point, link1, link2, frenet_trailer_error)

        else:
            # ==== OBSTACLES ====
            for obs in self.obstacles:
                obs_patch = plt.Circle((obs.x, obs.y), obs.radius, color='magenta', alpha=0.3)
                infl_patch = plt.Circle((obs.x, obs.y), obs.radius + self.clearance + self.influence_margin,
                                        color='magenta', fill=False, linestyle=':', alpha=0.5)
                ax_cartesian.add_patch(obs_patch)
                ax_cartesian.add_patch(infl_patch)

            mule_traj, = ax_cartesian.plot([], [], 'r-', label='Mule Path')
            trailer_traj, = ax_cartesian.plot([], [], color='lightblue', linestyle='-', label='Trailer Path')
            mule_point, = ax_cartesian.plot([], [], 'ro', label='Mule')
            hitch_point, = ax_cartesian.plot([], [], 'yo', label='Hitch')
            trailer_point, = ax_cartesian.plot([], [], 'ko', label='Trailer')
            link1, = ax_cartesian.plot([], [], 'r-', lw=1.5)
            link2, = ax_cartesian.plot([], [], 'k-', lw=1.5)
            ax_cartesian.legend()
            ax_cartesian.set_title("Cartesian Frame")

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
            for i in range(len(self.mules)):
                x, y = self.mules[i]
                x_ref, y_ref, theta_ref, kappa_ref = self.path_func(i * self.dt)
                # project mule on ref path: s increases linearly with v*dt
                s_accum += self.v * self.dt
                d = np.sin(theta_ref) * (x - x_ref) - np.cos(theta_ref) * (y - y_ref)
                frenet_coords.append([s_accum, d])
            frenet_coords = np.array(frenet_coords)

            def init():
                # Cartesian
                mule_traj.set_data([], [])
                trailer_traj.set_data([], [])
                mule_point.set_data([], [])
                hitch_point.set_data([], [])
                trailer_point.set_data([], [])
                link1.set_data([], [])
                link2.set_data([], [])
                # Frenet
                line_frenet.set_data([], [])
                mule_frenet_point.set_data([], [])
                return (mule_traj, mule_point, trailer_traj, hitch_point, trailer_point, link1, link2,
                        line_frenet, mule_frenet_point)

            def update(i):
                # === Cartesian ===
                mule = self.mules[i]
                hitch = self.hitches[i]
                trailer = self.trailers[i]

                mule_traj.set_data(self.mules[:i+1, 0], self.mules[:i+1, 1])
                trailer_traj.set_data(self.trailers[:i+1, 0], self.trailers[:i+1, 1])
                mule_point.set_data([mule[0]], [mule[1]])
                hitch_point.set_data([hitch[0]], [hitch[1]])
                trailer_point.set_data([trailer[0]], [trailer[1]])
                link1.set_data([mule[0], hitch[0]], [mule[1], hitch[1]])
                link2.set_data([hitch[0], trailer[0]], [hitch[1], trailer[1]])

                # === Frenet === (Only Mule)
                line_frenet.set_data(frenet_coords[:i+1, 0], frenet_coords[:i+1, 1])
                mule_frenet_point.set_data([frenet_coords[i, 0]], [frenet_coords[i, 1]])

                return (mule_traj, mule_point, trailer_traj, hitch_point, trailer_point, link1, link2,
                        line_frenet, mule_frenet_point)

        ani = animation.FuncAnimation(fig, update, frames=self.steps,
                                    init_func=init, interval=50, blit=True)
        # writer = animation.PillowWriter(fps=100) # Adjust fps (frames per second) as needed
        # ani.save('path_following.gif', writer=writer)
        # ani.save('obstacle_avoidance.gif', writer=writer)
        # ani.save('reverse_tracking.gif', writer=writer)
        self.displayTime(self.start_simulation, time.time())
        plt.show()

    def displayTime(self, start, end):
        k = int(end - start)
        mins = (k // 60)
        if end - start < 1:
            secs = (((end - start) * 10000) // 100) / 100
        else:
            secs = k - (mins * 60)
        print("Time taken: ", mins, "minutes, ", secs, "seconds")


if __name__ == "__main__":
    #========================================================================#
    #============= FORWARD FIGURE EIGHT WITH OBSTACLE AVOIDANCE =============#
    # direction toggle
    is_reverse = False

    path_func = lambda t: ParametricFunctions.figure_eight(t, a=10.0, v=1.0)
    controller = Controller(N=20, dt=0.1, v=1.0, is_reverse=is_reverse)
    sim = PathPlanningFrenetFrame(path_func, (0, 0, np.pi/4), controller, L=1.5, D=2.0, dt=0.1, T=80)

    # circular obstacle at (x,y) with radius r
    sim.add_obstacle(Obstacle(10.0, 0.0, 0.2),
                     clearance=1.0, influence_margin=1.0,
                     detour_span=5.0, pre_start=4.0, max_amplitude=0.5)
    sim.add_obstacle(Obstacle(-7.0, 5.0, 0.2),
                     clearance=1.0, influence_margin=1.0,
                     detour_span=5.0, pre_start=4.0, max_amplitude=0.5)

    sim.simulate()
    sim.animate()

    #========================================================================#
    #============ REVERSE FIGURE EIGHT WITHOUT OBSTACLE AVOIDANCE ===========#
    # direction toggle
    # is_reverse = True

    # path_func = lambda t: ParametricFunctions.figure_eight(t, a=10.0, v=1.0)
    # controller = Controller(N=20, dt=0.1, v=1.0, is_reverse=is_reverse)

    # x_start, y_start, theta_start, kappa_start = path_func(0.01)
    # start_pose = (x_start, y_start, CartesianFrenetConverter.normalize_angle(theta_start + np.pi))

    # sim = PathPlanningFrenetFrame(path_func, start_pose, controller, L=1.5, D=2.0, dt=0.1, T=80, is_reverse=is_reverse)

    # sim.simulate()
    # sim.animate()
