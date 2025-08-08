#!/Users/subhodeep/venv/bin/python

import math
import numpy as np
from scipy.optimize import minimize

class Controller:
    def __init__(self, N=12, dt=0.1, v=1.0,
                 w_d=50.0, w_heading=10.0, w_omega=0.01, w_omega_rate=1.0,
                 omega_bounds=(-2.0, 2.0), max_delta_omega=0.5, use_exp_smooth=False, smooth_alpha=0.3):
        """
        - N: MPC horizon
        - dt: timestep
        - v: forward speed (assumed constant)
        - weights: cost weights (tune these; high w_d forces tight lateral following)
        - omega_bounds: search bounds for omega
        - max_delta_omega: hard rate limit applied to final output (rad/s per control step)
        - use_exp_smooth: apply exponential smoothing on output (optional)
        """
        self.N = N
        self.dt = dt
        self.v = v

        # cost weights
        self.w_d = w_d
        self.w_heading = w_heading
        self.w_omega = w_omega
        self.w_omega_rate = w_omega_rate

        self.omega_min, self.omega_max = omega_bounds
        self.prev_omega = 0.0
        self.max_delta_omega = max_delta_omega

        # optional smoothing
        self.use_exp_smooth = use_exp_smooth
        self.smooth_alpha = smooth_alpha

    def step_dynamics(self, state, omega, L, D):
        x, y, theta, phi = state
        dx = self.v * np.cos(theta)
        dy = self.v * np.sin(theta)
        dtheta = omega
        dphi = (self.v * np.sin(theta - phi) - L * omega * np.cos(theta - phi)) / D
        return np.array([
            x + dx * self.dt,
            y + dy * self.dt,
            theta + dtheta * self.dt,
            phi + dphi * self.dt
        ])

    def frenet_errors(self, state, t, path_func):
        """
        Compute lateral deviation d and heading error (theta - theta_ref) based on a local path point.
        Uses the same simple projection used elsewhere in your code.
        """
        x, y, theta = state[0], state[1], state[2]
        rx, ry, theta_ref = path_func(t)
        # lateral deviation: sin(theta_ref)*(x - rx) - cos(theta_ref)*(y - ry)
        d = math.sin(theta_ref) * (x - rx) - math.cos(theta_ref) * (y - ry)
        # wrap heading difference to [-pi, pi]
        heading_err = math.atan2(math.sin(theta - theta_ref), math.cos(theta - theta_ref))
        return d, heading_err

    def cost(self, omega_seq, state0, t0, path_func, L, D):
        state = np.array(state0, dtype=float)
        J = 0.0
        prev_omega = self.prev_omega

        for i in range(self.N):
            omega = float(omega_seq[i])
            # forward simulate
            state = self.step_dynamics(state, omega, L, D)

            # time at this prediction step
            t_pred = t0 + (i + 1) * self.dt

            # compute Frenet errors (d, heading)
            d, heading_err = self.frenet_errors(state, t_pred, path_func)

            # accumulate cost:
            J += self.w_d * (d ** 2)
            J += self.w_heading * (heading_err ** 2)
            J += self.w_omega * (omega ** 2)
            J += self.w_omega_rate * ((omega - prev_omega) ** 2)

            prev_omega = omega

        return J

    def mpc(self, x, y, theta, phi, path_func, t0=0.0, L=1.5, D=2.0):
        state0 = [x, y, theta, phi]
        # initial guess: keep previous omega across the horizon for continuity
        init_guess = np.full(self.N, self.prev_omega, dtype=float)
        bounds = [(self.omega_min, self.omega_max)] * self.N

        res = minimize(self.cost, init_guess,
                       args=(state0, t0, path_func, L, D),
                       bounds=bounds, method='SLSQP',
                       options={'maxiter': 200, 'ftol': 1e-4})

        if res.success:
            omega_raw = float(res.x[0])
        else:
            # fallback: use previous omega
            omega_raw = self.prev_omega

        # optional exponential smoothing
        if self.use_exp_smooth:
            omega_out = self.smooth_alpha * omega_raw + (1.0 - self.smooth_alpha) * self.prev_omega
        else:
            omega_out = omega_raw

        # hard rate limit (per control step)
        maxd = self.max_delta_omega
        omega_out = float(np.clip(omega_out, self.prev_omega - maxd, self.prev_omega + maxd))

        # update previous
        self.prev_omega = omega_out
        return omega_out
