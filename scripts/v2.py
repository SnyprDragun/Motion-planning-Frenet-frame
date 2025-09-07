#!/Users/subhodeep/venv/bin/python
import time
import numpy as np
from tqdm import tqdm
from tabulate import tabulate
from dataclasses import dataclass
from typing import Tuple, Optional
import matplotlib.pyplot as plt
import matplotlib.animation as animation
try:
    from scipy.optimize import minimize
except ImportError as e:
    raise ImportError("This controller requires SciPy for optimization: pip install scipy") from e

'''
How does it work? Framework outline:

    First we have RobotDynamics:
        This class deals with the dynamics of the robot
        Parameters:
            - for number of trailers attached to mule (ideally recursive computation)
            - forward or reverse tracking
            - point body motion dynamics


    Then we have Path:
        This class is responsible for the final path to be traversed
        Parameters:
            - various parametric path functions
            - computes ideal path around obstcales


    Then we have FeedbackController:
        This class is the main control loop
        Parameters:
            - maintains a timer t in the background
            - checks current robot pose and ideal pose on path, and suggests control input on the dynamics at t
            - cost is only in frenet coordinate


    Finally we have PathPlanningFrenetFrame:
        This is our main class
        Parameters:
            - first chooses path and robot dynamics
            - tracks path for some t range
            - each step (t) consists of running FeedbackController on RobotDynamics to keep it on Path
            - plots trajectory in both cartesian and frenet



    Control LOOP
        - current state (all values known)
        - time step starts
        - first time step target pose for leading body
        - contoller return control action (velocity/acceleration) 
        - compute distance leading body covered in x-y and change in heading theta
            - control action isnt always optimal due to velocity bounds
        - reset robot pose from dynamics
            - depending on leading mass displacement, calculate pose of rest of robot
            - 
        - return leading body error to controller (in frenet)
        - repeat loop


        robot = Dynamics()
        for t in np.arange(0, T, timestep):
            lead_target = Path.circle(t)
            current_state = robot.current_pose
            control_action = Controller.mpc(current_state, lead_target)
            updated_lead_state = robot.update_lead_on_control(control_action)
            updated_robot_state = robot.update_robot_pose()
            store_all(current_state, control_action, updated_robot_state)

'''

class CartesianFrenetConverter:
    @staticmethod
    def cartesian_to_frenet():
        pass

    @staticmethod
    def frenet_to_cartesian():
        pass

    @ staticmethod
    def normalize_angle(angle):
        '''
        Normalize angle to [-pi, pi]
        '''
        return (angle + np.pi) % (2 * np.pi) - np.pi


class RobotDynamics:
    r'''
    Dynamics are given as follows:
        > The position of the robot (mule) is given by $(x, y, \theta)$. This point is at the centre of the front axle. 
        > Along the line from this point to the rear wheel and beyond is the hitch, at distance $L$.
        > Further, we have a trailer connected to this hitch point. 
        > The rear wheel of this trailer is at a distance $D$ from the hitch. 
        > We assume the line connecting the hitch and the two back wheels form an angle $\phi$.
    '''
    def __init__(self, offset, L_list, D_list, theta_list, phi_list, trailer_count=0, direction=True):
        self.L_list = L_list
        self.D_list = D_list
        self.theta_list = theta_list
        self.phi_list = phi_list

        self.trailer_count = trailer_count
        self.direction = direction

        self.mule_position = np.zeros(2)
        self.mule_orientation = self.theta_list[0]
        self.end_trailer_pose = offset

        self.control_actions = []
        self.calculate_mule_pose()

    def calculate_mule_pose(self):
        r'''
        If $(x_{trailer, N}, y_{trailer, N})$ is position of last trailer of the robot (offset), then mule position is calculated as:

            $x_{mule} = x_{trailer, N} + \sum_{i=1}^{N} D_{i}\cos(\phi_{i}) + L_{i}\cos(\theta_{i})$ and $y_{mule} = y_{trailer,N} + \sum_{i=1}^{N} D_{i}\sin(\phi_{i}) + L_{i}\sin(\theta_{i})$ where $N$ represents total numver of trailers
        '''
        self.mule_position = np.zeros(2)
        for i in range(self.trailer_count):
            xi = self.D_list[i] * np.cos(self.phi_list[i]) + self.L_list[i] * np.cos(self.theta_list[i])
            yi = self.D_list[i] * np.sin(self.phi_list[i]) + self.L_list[i] * np.sin(self.theta_list[i])
            self.mule_position = np.add(self.mule_position, np.array([xi, yi]))
        self.mule_position += self.end_trailer_pose
        return self.mule_position.tolist()

    def calculate_kth_hitch_trailer_pose(self, k):
        r'''
        If $(x_{mule}, y_{mule})$ is position of the mule, $k^{th}$ hitch and trailer position is calculated as:
            (assuming hitch-trailers are numbered $1$ to $N$, from mule to last trailer, and there are total $N$ hitch-trailer pairs)

            $x_{trailer, k} = x_{mule} - \sum_{i=1}^{k} D_{i}\cos(\phi_{i}) + L_{i}\cos(\theta_{i})$ and $y_{trailer, k} = y_{mule} - \sum_{i=1}^{k} D_{i}\sin(\phi_{i}) + L_{i}\sin(\theta_{i})$

            $x_{hitch, k} = x_{trailer, k} + D_{k}\cos(\phi_{k})$ and $y_{hitch, k} = y_{trailer, k} + D_{k}\sin(\phi_{k})$
        '''
        xi_trailer, yi_trailer = self.mule_position
        xi_hitch, yi_hitch = 0, 0

        try:
            for i in range(k):
                xi_trailer -= self.D_list[i] * np.cos(self.phi_list[i]) + self.L_list[i] * np.cos(self.theta_list[i])
                yi_trailer -= self.D_list[i] * np.sin(self.phi_list[i]) + self.L_list[i] * np.sin(self.theta_list[i])

            xi_hitch = xi_trailer + self.D_list[k - 1] * np.cos(self.phi_list[k - 1])
            yi_hitch = yi_trailer + self.D_list[k - 1] * np.sin(self.phi_list[k - 1])

            return [xi_hitch, yi_hitch], [xi_trailer, yi_trailer]

        except IndexError:
            print(f"ERROR: Robot has only {self.trailer_count} hitch-trailers, you are aksing pose for {k}")
            print("Proceeding with last hitch-trailer pose...")
            return self.calculate_kth_hitch_trailer_pose(self.trailer_count)

    def update_state(self, control_action, dt):
        r'''
        We approach this problem by first calculating for mule with $1$ hitch-trailer, then recursively calculate for the rest $N-1$ hitch-trailers.
        If control_action is $v$ and $\omega$, then we have $\dot x = v \cos(\theta)$, $\dot y = v \sin(\theta)$ and $\dot\theta = \omega$, where $(x, y, \theta)$ is mule pose ($\theta$ is angle at hitch).

        By kniematics, we have trailer dynamics as $v_{t} = v\cos(\theta-\phi) + L\omega\sin(\theta-\phi)$ and $\dot\phi = \frac{v\sin(\theta-\phi) - L\omega\cos(\theta-\phi)}{D}$

        So we can compute pose of mule and first hitch-trailer after each time step as:
        $x_{mule, new} = x_{mule, old} + v\cos(\theta) dt$
        $y_{mule, new} = y_{mule, old} + v\sin(\theta) dt$
        $\theta_{new} = \theta_{old} + \omega dt$
        $x_{trailer, new} = x_{trailer, old} + v_{t}\cos(\phi)dt$
        $y_{trailer, new} = y_{trailer, old} + v_{t}\sin(\phi)dt$
        $\phi_{new} = \phi_{old} + \dot\phi dt$

        For next hitch-trailer pair, current trailer acts as leading body, so their pose can be calculated in recursive fashion using the above formulae.
        Control action for new leading body is basically $v_{t}$ from previous step and $w_{next}$ which has to be calculated from previous $\phi$, similar to how $\dot\phi$ is calculated from $\theta$.
        '''
        v, omega = control_action
        x_mule_new = self.mule_position[0] + v * np.cos(self.theta_list[0]) * dt
        y_mule_new = self.mule_position[1] + v * np.sin(self.theta_list[0]) * dt

        new_theta_list = self.theta_list
        new_phi_list = self.phi_list
        v_ht , omega_ht = v, omega
        v_ht_new, omega_ht_new = 0, 0

        for i in range(self.trailer_count):
            self.control_actions.append([v_ht, omega_ht])
            _, trailer_old = self.calculate_kth_hitch_trailer_pose(i + 1)
            x_trailer_new = trailer_old[0] + self.v_trailer(v_ht, omega_ht, i) * np.cos(self.phi_list[i]) * dt
            y_trailer_new = trailer_old[1] + self.v_trailer(v_ht, omega_ht, i) * np.sin(self.phi_list[i]) * dt

            if i == self.trailer_count - 1:
                self.end_trailer_pose = [x_trailer_new, y_trailer_new]

            if i > 0:
                omega_ht_new = self.omega_next(v_ht, self.phi_dot(v_ht, omega_ht, i), i)
                v_ht_new = self.v_trailer(v_ht, omega_ht, i)
                omega_ht = omega_ht_new
                v_ht = v_ht_new

            new_theta_list[i] += omega_ht * dt
            new_phi_list[i] += self.phi_dot(v_ht, omega_ht, i) * dt

        self.theta_list = new_theta_list
        self.phi_list = new_phi_list
        self.mule_position = [x_mule_new, y_mule_new]
        self.mule_orientation = self.theta_list[0]

        return np.append(np.asarray(self.mule_position), np.asarray(self.mule_orientation))

    def v_trailer(self, v, omega, i):
        r'''
        $v_{t, i} = v\cos(\theta_{i}-\phi_{i}) + L_{i}\omega\sin(\theta_{i}-\phi_{i})$
        '''
        return v * np.cos(self.theta_list[i] - self.phi_list[i]) + self.L_list[i] * omega * np.sin(self.theta_list[i] - self.phi_list[i])

    def phi_dot(self, v, omega, i):
        r'''
        $\dot\phi_{t, i} = \frac{v\sin(\theta_{i}-\phi_{i}) - L_{i}\omega\cos(\theta_{i}-\phi_{i})}{D_{i}}$
        '''
        return (v * np.sin(self.theta_list[i] - self.phi_list[i]) - self.L_list[i] * omega * np.cos(self.theta_list[i] - self.phi_list[i])) / self.D_list[i]

    def omega_next(self, v_ht, phi_dot, i):
        r'''
        $\omega_{i} = \frac{v_{hitch, i-1} \sin(\phi_{i-1} - \theta_{i}) - D_{i-1} \dot\phi_{i-1} \cos(\phi_{i-1} - \theta_{i})}{L{i}}$
        '''
        return (v_ht * np.sin(self.phi_list[i-1] - self.theta_list[i]) - self.D_list[i-1] * phi_dot * np.cos(self.phi_list[i-1] - self.theta_list[i])) / self.L_list[i]

    def diagnostics(self):
        r'''
        Displays robot pose in tabular format
        +--------------+-------------+--------------+--------------+--------------+----....
        | Mule         | Hitch 1     | Trailer 1    | Hitch 2      | Trailer 2    |
        +==============+=============+==============+==============+==============+====....
        |[$x_{mule}$, $y_{mule}$] | [$x_{h, 1}$, $y_{h,1}$]  | [$x_{t, 1}$, $y_{t,1}$]   | [$x_{h, 2}$, $y_{h,2}$]   | [$x_{t, 2}$, $y_{t,2}$]   |
        +--------------+-------------+--------------+--------------+--------------+----....
        '''
        headers = ["Mule"]
        rows = [[[round(float(xy), 2) for xy in self.mule_position]]]
        for i in range(self.trailer_count):
            headers += [f"Hitch {i+1}"] + [f"Trailer {i+1}"]
            hitch, trailer = self.calculate_kth_hitch_trailer_pose(i+1)
            hitch = [round(float(xy), 2) for xy in hitch]
            trailer = [round(float(xy), 2) for xy in trailer]
            rows[0].append(hitch)
            rows[0].append(trailer)
        print(tabulate(rows, headers, tablefmt='grid'))


class Path():
    '''
    Parametric functions to generate path points. Also handles obstcale avoidance. 
    Returns. safe and optimal path point for traversal
    '''
    def __init__(self, shape):
        self.shape = shape
        self.obstacles = 0

    @staticmethod
    def circle(t, R=10.0, v=1.0):
        r'''
        Generates a circle path and its kinematic properties.
        Returns: $x$, $y$, $\theta$ (heading), $\kappa$ (curvature)
        • reaches back at starting point at $t = R * 2\pi$ (62.8)
        '''
        omega = v / R
        x = R * np.cos(omega * t)
        y = R * np.sin(omega * t)
        theta = np.arctan2(y, x) + np.pi / 2
        kappa = 1.0 / R
        return x, y, theta

    @staticmethod
    def ellipse(t, a=12.0, b=8.0, v=1.0):
        r'''
        Generates an ellipse path and its kinematic properties.
        Returns: $x$, $y$, $\theta$ (heading), $\kappa$ (curvature)
        '''
        omega = v / ((a + b) / 2.0)
        x = a * np.cos(omega * t)
        y = b * np.sin(omega * t)

        dx = -a * omega * np.sin(omega * t)
        dy = b * omega * np.cos(omega * t)
        ddx = -a * omega**2 * np.cos(omega * t)
        ddy = -b * omega**2 * np.sin(omega * t)

        theta = np.arctan2(dy, dx)
        speed_sq = dx**2 + dy**2

        if speed_sq < 1e-6:
            kappa = 0
        else:
            kappa = (dx * ddy - dy * ddx) / speed_sq**(1.5)
        return x, y, theta

    @staticmethod
    def figure_eight(t, a=10.0, v=1.0):
        r'''
        Generates a figure-eight path and its kinematic properties.
        Returns: $x$, $y$, $\theta$ (heading), $\kappa$ (curvature)
        '''
        omega = v / a

        dx = a * omega * np.cos(omega * t)
        dy = a * omega * np.cos(2 * omega * t)
        ddx = -a * omega**2 * np.sin(omega * t)
        ddy = -2 * a * omega**2 * np.sin(2 * omega * t)


        x = a * np.sin(omega * t)
        y = a * np.sin(omega * t) * np.cos(omega * t)
        theta = np.arctan2(dy, dx)
        speed_sq = dx**2 + dy**2

        if speed_sq < 1e-6:
            kappa = 0
        else:
            kappa = (dx * ddy - dy * ddx) / speed_sq**(1.5)
        return x, y, theta

    @staticmethod
    def straight_line(t, slope=0.0, intercept=0.0, v=1.0):
        r'''
        Generates a straight line path and its kinematic properties.
        Returns: $x$, $y$, $\theta$ (heading), $\kappa$ (curvature)
        '''
        x = v * t
        y = slope * x + intercept
        theta = np.arctan2(slope, 1.0)
        kappa = 0.0
        return x, y, theta

    def add_obstacle(self, *obstacles):
        '''
        Adds obstacle
        '''
        self.obstacles = obstacles
        pass

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


@dataclass
class Obstacle():
    x: float
    y: float
    radius: float


class PathPlanningFrenetFrame:
    def __init__(self, robot, target_path, controller, T, dt):
        self.robot = robot
        self.target_path = target_path
        self.controller = controller
        self.T = T
        self.dt = dt

        self.current_states = []
        self.target_states = []
        self.control_actions = []

        self.hitches = []
        self.trailers = []

        self.current_states_frenet = []
        self.target_states_frenet = []
        self.control_actions_frenet = []

    def store_hitch_trailer(self):
        hitches = []
        trailers = []

        for i in range(self.robot.trailer_count):
            hitch, trailer = self.robot.calculate_kth_hitch_trailer_pose(i + 1)
            hitches.append(hitch)
            trailers.append(trailer)

        self.hitches.append(hitches)
        self.trailers.append(trailers)

    def control_loop(self):
        current_state = np.append(self.robot.mule_position, self.robot.mule_orientation)

        for i in tqdm(range(int(self.T / self.dt)), desc="Simulating", unit="iterations"):
            t = i * self.dt

            target_state = Path.circle(t)
            control_action = self.controller.control(current_state, target_state)
            updated_state = robot.update_state(control_action, self.dt)
            current_state = updated_state

            self.store_hitch_trailer()
            self.current_states.append(current_state)
            self.target_states.append(target_state)
            self.control_actions.append(control_action)

    def diagnostics(self):
        r'''
        Displays control step in tabular format
        +---------------------+-----------------------+-------------------+---....
        |     Current Pose    |      Target Pose      |   Control Action  |
        +=====================+=======================+===================+===....
        | [$x_{mule}$, $y_{mule}$, $\theta_{mule}$] | [$x_{target}$, $y_{target}$, $\theta_{target}$] |       [$v$, $\omega$]      |
        +---------------------+-----------------------+-------------------+---....
        '''
        headers = ["Sl. no"] + ["Current Pose"] + ["Target Pose"] + ["Control Action"]
        rows = []
        for i in range(int(self.T / self.dt)):
            row = [
                i + 1,
                [f"{x:.4f}" for x in np.array(self.current_states[i]).flatten()],
                [f"{x:.4f}" for x in np.array(self.target_states[i]).flatten()],
                [f"{x:.4f}" for x in np.array(self.control_actions[i]).flatten()]
            ]
            rows.append(row)
        print(tabulate(rows, headers, tablefmt='grid'))

    def plot(self, interval=100):
        self.target_states = np.array(self.target_states)
        self.current_states = np.array(self.current_states)
        self.control_actions = np.array(self.control_actions)
        self.hitches = np.array(self.hitches)
        self.trailers = np.array(self.trailers)

        n_frames = int(self.T / self.dt)
        t = np.arange(n_frames)

        # --- Create layout with gridspec: big traj on left, two stacked on right
        fig = plt.figure(figsize=(12, 6))
        gs = fig.add_gridspec(2, 2, width_ratios=[2, 1])

        ax_traj = fig.add_subplot(gs[:, 0])   # trajectory big plot
        ax_v    = fig.add_subplot(gs[0, 1])  # v subplot
        ax_w    = fig.add_subplot(gs[1, 1])  # omega subplot

        # ---- Trajectory plot ----
        ax_traj.set_title("Trajectory")
        ax_traj.set_xlabel("X")
        ax_traj.set_ylabel("Y")
        ax_traj.plot(self.target_states[:, 0], self.target_states[:, 1], "k--", label="Target")
        current_line, = ax_traj.plot([], [], "r-", label="Current")
        ax_traj.legend()
        ax_traj.set_aspect("equal", adjustable="datalim")

        # ---- Control plots ----
        ax_v.set_title("Linear Velocity v")
        ax_v.set_xlabel("Time")
        ax_v.set_ylabel("v")
        line_v, = ax_v.plot([], [], "b-")

        ax_w.set_title("Angular Velocity ω")
        ax_w.set_xlabel("Time")
        ax_w.set_ylabel("ω")
        line_w, = ax_w.plot([], [], "g-")

        # ---- Hitch & Trailer plots ----
        n_pairs = self.hitches.shape[1]
        mule_scats = [ax_traj.plot([], [], "ro", label=f"Mule")[0]]
        hitch_scats = [ax_traj.plot([], [], "bo", label=f"Hitch {i+1}" if i == 0 else "")[0]
                    for i in range(n_pairs)]
        trailer_scats = [ax_traj.plot([], [], "go", label=f"Trailer {i+1}" if i == 0 else "")[0]
                        for i in range(n_pairs)]

        # --- Update function ---
        def update(frame):
            # trajectory + control
            current_line.set_data(self.current_states[:frame, 0], self.current_states[:frame, 1])
            line_v.set_data(t[:frame], self.control_actions[:frame, 0])
            line_w.set_data(t[:frame], self.control_actions[:frame, 1])

            # hitches & trailers
            for i in range(n_pairs):
                mx, my = self.current_states[frame-1:frame, 0], self.current_states[frame-1:frame, 1]
                hx, hy = self.hitches[frame, i]
                tx, ty = self.trailers[frame, i]
                mule_scats[0].set_data([mx], [my])
                hitch_scats[i].set_data([hx], [hy])
                trailer_scats[i].set_data([tx], [ty])

            # rescale
            ax_traj.relim(); ax_traj.autoscale_view()
            ax_v.relim(); ax_v.autoscale_view()
            ax_w.relim(); ax_w.autoscale_view()

            return [current_line, line_v, line_w] + mule_scats + hitch_scats + trailer_scats

        ani = animation.FuncAnimation(fig, update, frames=n_frames, interval=interval, blit=False, repeat=False)
        plt.show()

    def display_time(self, start, end):
        k = int(end - start)
        mins = (k // 60)
        if end - start < 1:
            secs = (((end - start) * 10000) // 100) / 100
        else:
            secs = k - (mins * 60)
        print(f"Time taken: {mins} minutes {secs} seconds")


if __name__ == "__main__":
    r'''
    < ================= Main ================= >

    Choose:
        > robot dynamics - no. of trailers, direction of motion
        > controller - mpc or other
        > path - straight line, circle, ellipse, figure-eight

    Specify:
        > initial position of last trailer (offset)
        > initial values of all angles ($\theta$, $\phi$)
        > lengths of joining rods ($L$, $D$)
    '''

    start = time.time()
    robot = RobotDynamics([10,-10.5], [1.5, 1.5, 1.5], [2.0, 2.0, 2.0], [np.pi/2, np.pi/2, np.pi/2], [np.pi/2, np.pi/2, np.pi/2], trailer_count=3, direction=True)
    # robot = RobotDynamics([10,-7], [1.5, 1.5], [2.0, 2.0], [np.pi/2, np.pi/2], [np.pi/2, np.pi/2], trailer_count=2, direction=True)
    # robot = RobotDynamics([-3.5*np.pi/4, -3.5*np.pi/4], [1.5], [2.0], [np.pi/4], [np.pi/4], trailer_count=1, direction=True)
    # robot = RobotDynamics([10,-3.5], [1.5], [2.0], [np.pi/2], [np.pi/2], trailer_count=1, direction=True)
    robot.diagnostics()
    mpc = UnicycleMPC(MPCParams(dt=0.1, N=20))
    circle = Path("circle")
    # circle.add_obstacles()
    trajectory = PathPlanningFrenetFrame(robot=robot, target_path=circle, controller=mpc, T=65, dt=0.1)
    trajectory.control_loop()
    # trajectory.diagnostics()
    robot.diagnostics()
    trajectory.display_time(start, time.time())
    trajectory.plot()
