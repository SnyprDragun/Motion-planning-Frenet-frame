#!/Users/subhodeep/venv/bin/python

import math
import numpy as np
from scipy.optimize import minimize

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
