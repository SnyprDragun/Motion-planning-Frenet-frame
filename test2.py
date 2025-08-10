#!/Users/subhodeep/venv/bin/python

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.optimize import minimize

class CartesianFrenetConverter:
    """
    A class for converting states between Cartesian and Frenet coordinate systems.
    This class remains unchanged as its functionality is universal.
    """

    @staticmethod
    def cartesian_to_frenet(rs, rx, ry, rtheta, rkappa, rdkappa, x, y, v, a, theta, kappa):
        """
        Convert state from Cartesian coordinate to Frenet coordinate
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
        # Ensure one_minus_kappa_r_d is not zero to avoid division by zero
        if abs(one_minus_kappa_r_d) < 1e-6:
            return None, None

        d_prime = one_minus_kappa_r_d * tan_delta_theta

        kappa_r_d_prime = rdkappa * d + rkappa * d_prime

        d_pprime = (-kappa_r_d_prime * tan_delta_theta +
                  one_minus_kappa_r_d / (cos_delta_theta * cos_delta_theta) *
                  (kappa * one_minus_kappa_r_d / cos_delta_theta - rkappa))

        s = rs
        s_dot = v * cos_delta_theta / one_minus_kappa_r_d

        delta_theta_prime = one_minus_kappa_r_d / cos_delta_theta * kappa - rkappa
        s_ddot = (a * cos_delta_theta -
                  s_dot * s_dot *
                  (d_prime * delta_theta_prime - kappa_r_d_prime)) / one_minus_kappa_r_d

        return [s, s_dot, s_ddot], [d, d_prime, d_pprime]

    @staticmethod
    def normalize_angle(angle):
        """
        Normalize angle to [-pi, pi]
        """
        a = math.fmod(angle + math.pi, 2.0 * math.pi)
        if a < 0.0:
            a += 2.0 * math.pi
        return a - math.pi


class ParametricFunctions:
    """
    Generates parametric paths. Now also returns the path's curvature.
    """
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

class PathPlanningFrenetFrame:
    """
    Main simulation class for reverse path tracking.
    """
    def __init__(self, path_func, start_trailer_pose, controller, L=1.5, D=2.0, dt=0.1, T=60):
        self.path_func = path_func
        self.controller = controller
        self.L = L
        self.D = D
        self.dt = dt
        self.T = T
        self.steps = int(T / dt)
        
        self.x, self.y, self.theta = start_trailer_pose
        self.phi = 0.0
        self.v = controller.v

        # Data logging lists
        self.mule_path = []
        self.hitch_path = []
        self.trailer_path = []
        self.frenet_trailer_path = []
        
        # Pre-compute a high-resolution reference path for accurate plotting
        self.ref_path_dense = None
        self.ref_s_dense = None
        self._precompute_ref_path()

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

    def simulate(self):
        """Run the full reverse tracking simulation."""
        last_closest_idx = 0
        for i in range(self.steps):
            t = i * self.dt
            
            omega = self.controller.mpc(self.x, self.y, self.theta, self.phi,
                                        self.path_func, t0=t, L=self.L, D=self.D)

            # Update state using reverse dynamics
            theta_m = CartesianFrenetConverter.normalize_angle(self.theta + self.phi)
            vx_h = -self.v * math.cos(theta_m) + self.L * omega * math.sin(theta_m)
            vy_h = -self.v * math.sin(theta_m) - self.L * omega * math.cos(theta_m)
            v_t = vx_h * math.cos(self.theta) + vy_h * math.sin(self.theta)
            v_perp = -vx_h * math.sin(self.theta) + vy_h * math.cos(self.theta)
            dtheta_dt = v_perp / self.D
            dphi_dt = omega - dtheta_dt
            
            self.x += v_t * math.cos(self.theta) * self.dt
            self.y += v_t * math.sin(self.theta) * self.dt
            self.theta = CartesianFrenetConverter.normalize_angle(self.theta + dtheta_dt * self.dt)
            self.phi = CartesianFrenetConverter.normalize_angle(self.phi + dphi_dt * self.dt)

            # Log Cartesian poses for plotting
            mule, hitch = self.compute_mule_hitch_from_trailer()
            self.trailer_path.append([self.x, self.y])
            self.hitch_path.append(hitch)
            self.mule_path.append(mule)

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
            d = math.copysign(math.hypot(dx, dy), math.sin(rtheta) * dx - math.cos(rtheta) * dy)
            
            self.frenet_trailer_path.append([s_true, d])


        # Convert lists to numpy arrays
        self.trailer_path = np.array(self.trailer_path)
        self.hitch_path = np.array(self.hitch_path)
        self.mule_path = np.array(self.mule_path)
        self.frenet_trailer_path = np.array(self.frenet_trailer_path)

    def animate(self):
        """Animate the simulation results with an added error plot."""
        fig = plt.figure(figsize=(16, 8))
        ax_cartesian = fig.add_subplot(1, 2, 1)
        ax_frenet = fig.add_subplot(1, 2, 2)

        # === Cartesian Plot ===
        ax_cartesian.set_aspect('equal')
        ax_cartesian.set_xlim(-15, 15)
        ax_cartesian.set_ylim(-15, 15)
        ax_cartesian.plot(self.ref_path_dense[:, 0], self.ref_path_dense[:, 1], 'gray', linestyle='--', label='Reference Path')
        line_trailer, = ax_cartesian.plot([], [], 'b-', label='Trailer Path')
        trailer_point, = ax_cartesian.plot([], [], 'bo', ms=8, label='Trailer (Lead)')
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
        line_frenet_error, = ax_frenet.plot([], [], 'g-', label='Lateral Error (d)')
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
            line_frenet_error.set_data([], [])
            return (line_trailer, trailer_point, hitch_point, mule_point, link1, link2, line_frenet_error)

        def update(i):
            trailer = self.trailer_path[i]
            hitch = self.hitch_path[i]
            mule = self.mule_path[i]

            line_trailer.set_data(self.trailer_path[:i+1, 0], self.trailer_path[:i+1, 1])
            trailer_point.set_data([trailer[0]], [trailer[1]])
            hitch_point.set_data([hitch[0]], [hitch[1]])
            mule_point.set_data([mule[0]], [mule[1]])
            link1.set_data([trailer[0], hitch[0]], [trailer[1], hitch[1]])
            link2.set_data([hitch[0], mule[0]], [hitch[1], mule[1]])

            # Corrected Frenet plot update
            if len(self.frenet_trailer_path) > i:
                frenet_data = self.frenet_trailer_path[:i+1]
                line_frenet_error.set_data(frenet_data[:, 0], frenet_data[:, 1])
            
            return (line_trailer, trailer_point, hitch_point, mule_point, link1, link2, line_frenet_error)

        ani = animation.FuncAnimation(fig, update, frames=self.steps,
                                    init_func=init, interval=50, blit=True)
        plt.show()


if __name__ == "__main__":
    # The parametric function now defines the desired path for the TRAILER
    path_func = lambda t: ParametricFunctions.figure_eight(t, a=10.0, v=1.0)

    # Controller with fine-tuned weights for smooth, stable motion
    controller = Controller(N=20, dt=0.1, v=1.0, is_reverse=True)
    
    # Get the initial pose from the path function at t=0.01 to avoid singularity
    x_start, y_start, theta_start, _ = path_func(0.01)
    # For reverse tracking, the trailer must face the opposite direction of the path's tangent
    start_pose_trailer = (x_start, y_start, CartesianFrenetConverter.normalize_angle(theta_start + math.pi))
    
    # Initialize and run the simulation
    sim = PathPlanningFrenetFrame(path_func, start_pose_trailer, controller, L=1.5, D=2.0, T=80)
    print("Starting reverse tracking simulation with improved controller...")
    sim.simulate()
    print("Simulation complete. Starting animation...")
    sim.animate()
