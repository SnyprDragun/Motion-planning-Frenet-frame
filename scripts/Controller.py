#!/Users/subhodeep/venv/bin/python

import numpy as np
from scipy.optimize import minimize
from CartesianFrenetConverter import CartesianFrenetConverter


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
