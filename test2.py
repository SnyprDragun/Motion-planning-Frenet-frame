#!/Users/subhodeep/venv/bin/python
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.optimize import minimize
from dataclasses import dataclass

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
            raise ValueError("The reference point s and s_condition[0] don't match")

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


class Controller:
    def __init__(self, N=12, dt=0.1, v=1.0,
                 w_d_mule=1200.0, w_d_hitch=0.0, w_d_trailer=0.0,
                 w_omega=0.01, w_omega_rate=1.0,
                 omega_bounds=(-2.0, 2.0), max_delta_omega=0.5,
                 use_exp_smooth=False, smooth_alpha=0.3):
        """
        MPC controller using only lateral Frenet errors of mule, hitch, trailer.
        Key tuning names:
        - w_d_mule: lateral deviation weight for mule (increase to follow path more strictly)
        - controller lookahead = N * dt
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


@dataclass
class Obstacle:
    x: float
    y: float
    radius: float


class PathPlanningFrenetFrame:
    """
    Frenet-based detour using a Gaussian (bell-curve) bump d(s) precomputed and latched
    onto the parametric path.
    """
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

        # Precompute Frenet reference path logs
        self.ref_s = []
        self.ref_d = []

        # Obstacles list
        self.obstacles = []
        # clearance around obstacle surface (how far laterally we want to be)
        self.clearance = 1.0
        # influence margin beyond (radius + clearance) to trigger detection
        self.influence_margin = 3.0
        # detour span (s-length) over which Gaussian bump is significant
        self.detour_span = 4.0  # meters; tune if needed
        # Gaussian sigma will be set as detour_span / 3 (so tails near zero at ends)

        # Frenet detour state
        self.detour_active = False
        self.detour_s_center = None
        self.detour_s0 = None
        self.detour_s1 = None
        self.detour_amplitude = None
        self.detour_sigma = None
        self.detour_obs = None

        # debug flag
        self.debug = False

        # sim time
        self.sim_time = 0.0

    def add_obstacle(self, obstacle: Obstacle, clearance: float = 1.0, influence_margin: float = 3.0, detour_span: float = 4.0):
        self.obstacles.append(obstacle)
        self.clearance = clearance
        self.influence_margin = influence_margin
        self.detour_span = detour_span

    def compute_hitch_trailer(self):
        x_h = self.x - self.L * np.cos(self.theta)
        y_h = self.y - self.L * np.sin(self.theta)
        x_t = x_h - self.D * np.cos(self.phi)
        y_t = y_h - self.D * np.sin(self.phi)
        return np.array([x_h, y_h]), np.array([x_t, y_t])

    def _s_of_time(self, t):
        """Arc-length along parametric path used as reference s (we used s = v * t)."""
        return t * self.v

    def _path_point(self, t):
        """Return parametric path point (x,y,theta) at time t."""
        return self.path_func(t)

    def _distance_point_to_obs(self, x, y, obs):
        return math.hypot(x - obs.x, y - obs.y)

    def _should_plan_detour(self, current_t):
        """
        Scan the parametric path forward over MPC lookahead; if any point is within
        obs.radius + clearance + influence_margin, request detour planning.
        """
        lookahead = self.controller.N * self.controller.dt
        steps = max(1, int(lookahead / self.dt))
        for k in range(steps):
            t_check = current_t + k * self.dt
            px, py, _ = self._path_point(t_check)
            for obs in self.obstacles:
                inner = obs.radius + self.clearance
                outer = inner + self.influence_margin
                d = self._distance_point_to_obs(px, py, obs)
                if d <= outer:
                    if self.debug:
                        print(f"[detect] at t={current_t:.2f} horizon found t_check={t_check:.2f} d={d:.2f} obs@({obs.x},{obs.y})")
                    return obs, t_check
        return None, None

    def _find_s_center(self, detect_t, obs, search_span=6.0):
        """
        Find t in [detect_t - search_span, detect_t + search_span] that minimizes distance to obstacle.
        Use coarse sampling to find s_center.
        """
        t0 = max(0.0, detect_t - search_span)
        t1 = detect_t + search_span
        N = max(20, int((t1 - t0) / self.dt))
        best_t = detect_t
        best_d = float('inf')
        for i in range(N + 1):
            tt = t0 + (t1 - t0) * (i / N)
            px, py, _ = self._path_point(tt)
            d = self._distance_point_to_obs(px, py, obs)
            if d < best_d:
                best_d = d
                best_t = tt
        s_center = self._s_of_time(best_t)
        if self.debug:
            print(f"[find_center] best_t={best_t:.2f} best_d={best_d:.2f} s_center={s_center:.2f}")
        return s_center, best_d

    def _plan_detour_gaussian(self, detect_t, obs):
        """
        Plan a Gaussian detour bump in d(s):
            d(s) = A * exp(- (s - s_center)^2 / (2 * sigma^2) )
        We set:
            - s_center at closest parametric point to obstacle
            - span = detour_span, so set sigma = detour_span / 3 (≈ 99% mass inside ~±1.5 sigma)
            - amplitude A chosen to move away from obstacle: A = sign * (obs.radius + clearance)
        """
        s_center, best_d = self._find_s_center(detect_t, obs)

        half_span = 0.5 * self.detour_span
        s0 = max(0.0, s_center - half_span)
        s1 = s_center + half_span

        # Determine sign: compute normal left-vector at center and dot with obstacle vector
        t_center = s_center / self.v
        rx, ry, theta_ref = self._path_point(t_center)
        vx = obs.x - rx
        vy = obs.y - ry
        # left normal (points to left of path)
        nx = -math.sin(theta_ref)
        ny = math.cos(theta_ref)
        dot = vx * nx + vy * ny
        # if obstacle is left (dot>0), detour should be to right (negative d)
        sign = -1.0 if dot > 0.0 else 1.0
        A = sign * (obs.radius + self.clearance)

        sigma = max(1e-3, self.detour_span / 3.0)

        # Store detour parameters
        self.detour_active = True
        self.detour_s_center = s_center
        self.detour_s0 = s0
        self.detour_s1 = s1
        self.detour_amplitude = A
        self.detour_sigma = sigma
        self.detour_obs = obs

        if self.debug:
            print(f"[plan_gauss] s0={s0:.2f}, s_center={s_center:.2f}, s1={s1:.2f}, A={A:.2f}, sigma={sigma:.2f}")

    def _maybe_deactivate_detour(self, current_s):
        """Deactivate detour once we have passed s1."""
        if not self.detour_active:
            return
        if current_s > self.detour_s1 + 1e-2:
            if self.debug:
                print(f"[detour_done] current_s={current_s:.2f} > s1={self.detour_s1:.2f}, deactivating detour")
            self.detour_active = False
            self.detour_s_center = None
            self.detour_s0 = None
            self.detour_s1 = None
            self.detour_amplitude = None
            self.detour_sigma = None
            self.detour_obs = None

    def _gaussian_bump(self, s):
        """
        Evaluate Gaussian bump and its first two derivatives wrt s:
            d(s) = A * exp(- (s - mu)^2 / (2 sigma^2))
            d'(s) = A * (-(s-mu) / sigma^2) * exp(...)
            d''(s) = A * ( ((s-mu)^2 / sigma^4) - (1 / sigma^2) ) * exp(...)
        Returns (d, d_s, d_ss)
        """
        if (not self.detour_active) or (self.detour_s_center is None):
            return 0.0, 0.0, 0.0
        mu = self.detour_s_center
        sigma = self.detour_sigma
        A = self.detour_amplitude
        z = (s - mu)
        expv = math.exp(-0.5 * (z * z) / (sigma * sigma))
        d = A * expv
        d_s = A * (-z / (sigma * sigma)) * expv
        d_ss = A * ((z * z) / (sigma * sigma * sigma * sigma) - 1.0 / (sigma * sigma)) * expv
        return d, d_s, d_ss

    def avoidance_path_frenet(self, t_pred):
        """
        If detour is active and s_pred near detour, return frenet->cartesian conversion
        using the Gaussian lateral offset. Otherwise return parametric point.
        """
        rx, ry, rtheta = self._path_point(t_pred)
        s = self._s_of_time(t_pred)

        # If no detour, return parametric
        if not self.detour_active:
            return rx, ry, rtheta

        # If outside planned detour support, return parametric
        if s < (self.detour_s0 - 1e-8) or s > (self.detour_s1 + 1e-8):
            return rx, ry, rtheta

        # Compute d(s), d'(s), d''(s)
        d, d_s, d_ss = self._gaussian_bump(s)

        # s_condition: [s, s_dot, s_ddot] ; s_dot = v ; s_ddot = 0
        s_cond = [s, self.v, 0.0]
        d_cond = [d, d_s, d_ss]

        # Convert frenet to cartesian using reference at same s (t_pred)
        try:
            x, y, theta_det, kappa, v_out, a_out = CartesianFrenetConverter.frenet_to_cartesian(
                rs=s, rx=rx, ry=ry, rtheta=rtheta, rkappa=0.0, rdkappa=0.0,
                s_condition=s_cond, d_condition=d_cond)
        except Exception as e:
            if self.debug:
                print(f"[frenet_conv_err] {e}")
            return rx, ry, rtheta
        return x, y, theta_det

    def simulate(self):
        for i in range(self.steps):
            t = i * self.dt
            current_s = self._s_of_time(t)

            # If no detour active, check if we need to plan one
            if not self.detour_active:
                obs, detect_t = self._should_plan_detour(t)
                if obs is not None and detect_t is not None:
                    self._plan_detour_gaussian(detect_t, obs)

            # Deactivate if we passed the detour
            self._maybe_deactivate_detour(current_s)

            # Build path_func closure for MPC (uses frenet detour if active)
            path_func_local = lambda t_pred: self.avoidance_path_frenet(t_pred)

            omega = self.controller.mpc(self.x, self.y, self.theta, self.phi,
                                        path_func_local, t0=t, L=self.L, D=self.D)

            # dynamics
            self.phi += ((self.v * np.sin(self.theta - self.phi) - self.L * omega * np.cos(self.theta - self.phi)) / self.D) * self.dt
            self.x += self.v * np.cos(self.theta) * self.dt
            self.y += self.v * np.sin(self.theta) * self.dt
            self.theta += omega * self.dt
            self.theta = CartesianFrenetConverter.normalize_angle(self.theta)

            # advance sim time
            self.sim_time += self.dt

            hitch, trailer = self.compute_hitch_trailer()

            # Log Cartesian
            self.states.append([self.x, self.y])
            self.hitches.append(hitch)
            self.trailers.append(trailer)

            # Log Frenet (for mule, hitch, trailer) using original parametric reference
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
            self.ref_s.append(current_s)
            self.ref_d.append(0.0)

        # convert lists to arrays
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

        # plot obstacles
        for obs in self.obstacles:
            obs_patch = plt.Circle((obs.x, obs.y), obs.radius, color='magenta', alpha=0.3)
            infl_patch = plt.Circle((obs.x, obs.y), obs.radius + self.clearance + self.influence_margin,
                                    color='magenta', fill=False, linestyle=':', alpha=0.5)
            ax_cart.add_patch(obs_patch)
            ax_cart.add_patch(infl_patch)

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
            s_accum += self.v * self.dt
            d = np.sin(theta_ref) * (x - x_ref) - np.cos(theta_ref) * (y - y_ref)
            frenet_coords.append([s_accum, d])
        frenet_coords = np.array(frenet_coords)

        def init():
            line_traj.set_data([], [])
            mule_point.set_data([], [])
            hitch_point.set_data([], [])
            trailer_point.set_data([], [])
            link1.set_data([], [])
            link2.set_data([], [])
            line_frenet.set_data([], [])
            mule_frenet_point.set_data([], [])
            return (line_traj, mule_point, hitch_point, trailer_point, link1, link2,
                    line_frenet, mule_frenet_point)

        def update(i):
            mule = self.states[i]
            hitch = self.hitches[i]
            trailer = self.trailers[i]

            line_traj.set_data(self.states[:i+1, 0], self.states[:i+1, 1])
            mule_point.set_data([mule[0]], [mule[1]])
            hitch_point.set_data([hitch[0]], [hitch[1]])
            trailer_point.set_data([trailer[0]], [trailer[1]])
            link1.set_data([mule[0], hitch[0]], [mule[1], hitch[1]])
            link2.set_data([hitch[0], trailer[0]], [hitch[1], trailer[1]])

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
    sim = PathPlanningFrenetFrame(path_func, (0, 0, np.pi/4), controller, L=1.5, D=2.0, dt=0.1, T=60)

    # Add a circular obstacle at (10,0) radius 1
    sim.add_obstacle(Obstacle(10.0, 0.0, 0.5), clearance=0.5, influence_margin=0.5, detour_span=4.0)
    sim.add_obstacle(Obstacle(-10.0, 0.0, 0.5), clearance=0.5, influence_margin=0.5, detour_span=4.0)

    # Optional debugging
    # sim.debug = True
    # controller.w_d_mule = 2000.0  # tighten lateral following if you want

    sim.simulate()
    sim.animate()
