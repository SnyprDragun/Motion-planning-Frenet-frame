#!/Users/subhodeep/venv/bin/python3
'''classes for MPC'''
import numpy as np
from dataclasses import dataclass
from typing import Tuple, Optional
from CartesianFrenetConverter import CartesianFrenetConverter
try:
    from scipy.optimize import minimize
except ImportError as e:
    raise ImportError("This controller requires SciPy for optimization: pip install scipy") from e

@dataclass
class MPCParams:
    dt: float = 0.1                  # sampling time [s]
    N: int = 15                      # prediction horizon
    v_min: float = -0.8              # min linear velocity [m/s]
    v_max: float = 1.0               # max linear velocity [m/s]
    w_min: float = -2.5              # min angular velocity [rad/s]
    w_max: float = 2.5               # max angular velocity [rad/s]
    dv_max: float = 0.5              # max |Δv| per step [m/s]
    dw_max: float = 1.0              # max |Δω| per step [rad/s]
    # costs
    q_xy: float = 2.0                # weight on (x,y) position error
    q_th: float = 0.3                # weight on heading error
    r_v: float = 0.02                # weight on absolute v
    r_w: float = 0.02                # weight on absolute ω
    s_v: float = 0.5                 # weight on |Δv|
    s_w: float = 0.3                 # weight on |Δω|
    qT_xy: float = 6.0               # terminal (x,y)
    qT_th: float = 1.0               # terminal heading


class UnicycleMPC:
    """
    Model Predictive Controller for a unicycle robot.
    State:   x = [x, y, theta]
    Control: u = [v, w]  (linear, angular)
    Dynamics:
        x_{k+1} = x_k + dt * [ v*cos(theta_k), v*sin(theta_k), w ]
    """
    def __init__(self, params: Optional[MPCParams] = None):
        self.p = params or MPCParams()
        # Keep last solution for warm-start
        self._u_prev_seq = np.zeros(2 * self.p.N)  # [v0..vN-1, w0..wN-1]

    # ----------------- public API -----------------
    def control(self, x_now: np.ndarray, x_target: np.ndarray) -> Tuple[float, float]:
        """
        Compute the first control action (v, w).
        x_now:     [x, y, theta]
        x_target:  [x*, y*, theta*]  (theta* optional: if NaN, ignore heading error)
        """
        x_now = np.asarray(x_now, dtype=float).reshape(3)
        x_ref = np.asarray(x_target, dtype=float).reshape(3)
        theta_ref_valid = np.isfinite(x_ref[2])

        # Warm start: shift previous optimal controls forward one step
        u0 = self._warm_start()

        # Box bounds on controls
        bounds = [(self.p.v_min, self.p.v_max)] * self.p.N + [(self.p.w_min, self.p.w_max)] * self.p.N

        # Build inequality constraints for slew-rate (Δu) limits: |u_k - u_{k-1}| <= d_max
        A_ineq, b_ineq = self._delta_constraints()

        cons = []
        if A_ineq is not None:
            cons.append({
                "type": "ineq",
                "fun": lambda u, A=A_ineq, b=b_ineq: b - A @ u
            })

        # Objective
        def objective(u: np.ndarray) -> float:
            v_seq = u[:self.p.N]
            w_seq = u[self.p.N:]
            # simulate
            traj = self._rollout(x_now, v_seq, w_seq)
            # stage costs
            pos_err = traj[:, :2] - x_ref[:2]          # (N+1, 2)
            e_xy2 = np.sum(pos_err**2, axis=1)        # (N+1,)
            if theta_ref_valid:
                e_th = np.array([CartesianFrenetConverter.normalize_angle(th - x_ref[2]) for th in traj[:, 2]])
                e_th2 = e_th**2
            else:
                e_th2 = np.zeros(self.p.N + 1)

            # smoothness (Δu) and effort
            dv = np.diff(v_seq, prepend=v_seq[0])
            dw = np.diff(w_seq, prepend=w_seq[0])

            J_stage = np.sum(self.p.q_xy * e_xy2[:-1] + self.p.q_th * e_th2[:-1]
                             + self.p.r_v * (v_seq**2) + self.p.r_w * (w_seq**2)
                             + self.p.s_v * (np.abs(dv)) + self.p.s_w * (np.abs(dw)))

            # terminal
            J_term = self.p.qT_xy * e_xy2[-1] + self.p.qT_th * e_th2[-1]
            return J_stage + J_term

        res = minimize(
            objective, u0, method="SLSQP", bounds=bounds, constraints=cons,
            options={"maxiter": 100, "ftol": 1e-4, "disp": False}
        )

        u_opt = res.x if res.success else u0
        self._u_prev_seq = u_opt.copy()

        v_cmd = float(np.clip(u_opt[0], self.p.v_min, self.p.v_max))
        w_cmd = float(np.clip(u_opt[self.p.N], self.p.w_min, self.p.w_max))
        return [v_cmd, w_cmd]

    # ----------------- helpers -----------------
    def _rollout(self, x0: np.ndarray, v_seq: np.ndarray, w_seq: np.ndarray) -> np.ndarray:
        """Simulate N steps; return (N+1, 3) states including x0."""
        dt = self.p.dt
        x = x0.copy()
        traj = [x.copy()]
        for k in range(self.p.N):
            v = v_seq[k]
            w = w_seq[k]
            x[0] = x[0] + dt * v * np.cos(x[2])
            x[1] = x[1] + dt * v * np.sin(x[2])
            x[2] = CartesianFrenetConverter.normalize_angle(x[2] + dt * w)
            traj.append(x.copy())
        return np.vstack(traj)

    def _warm_start(self) -> np.ndarray:
        """Shift previous sequence forward; keep last value for the tail."""
        u = self._u_prev_seq
        N = self.p.N
        v_seq = u[:N]
        w_seq = u[N:]
        v_ws = np.r_[v_seq[1:], v_seq[-1]]
        w_ws = np.r_[w_seq[1:], w_seq[-1]]
        return np.r_[v_ws, w_ws]

    def _delta_constraints(self):
        """
        Build |Δv_k|<=dv_max, |Δw_k|<=dw_max as linear inequalities:
            -dv_max <= v_k - v_{k-1} <= dv_max   (similar for w)
        SLSQP accepts A u <= b. We express both sides by stacking.
        Use first element relative to previous optimal first element (fixed at itself -> no constraint tightening).
        """
        N = self.p.N
        dv = self.p.dv_max
        dw = self.p.dw_max
        if not np.isfinite(dv) and not np.isfinite(dw):
            return None, None

        # Variables: u = [v0..vN-1, w0..wN-1]
        dim = 2 * N
        rows = []
        rhs = []

        # v deltas
        if np.isfinite(dv):
            # k=0: we don't constrain against unknown previous command; instead cap |v0| change by bounds already.
            for k in range(1, N):
                row = np.zeros(dim)
                row[k] = 1.0      # v_k
                row[k - 1] = -1.0 # -v_{k-1}
                rows.append(row.copy());  rhs.append(dv)   #  v_k - v_{k-1} <= dv
                rows.append(-row.copy()); rhs.append(dv)   # -v_k + v_{k-1} <= dv  => -(v_k - v_{k-1}) <= dv

        # w deltas
        if np.isfinite(dw):
            offset = N
            for k in range(1, N):
                row = np.zeros(dim)
                row[offset + k] = 1.0
                row[offset + k - 1] = -1.0
                rows.append(row.copy());  rhs.append(dw)
                rows.append(-row.copy()); rhs.append(dw)

        if not rows:
            return None, None
        A = np.vstack(rows)
        b = np.asarray(rhs)
        return A, b
