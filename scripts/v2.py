#!/Users/subhodeep/venv/bin/python

import numpy as np
from tabulate import tabulate
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
    pass

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
        self.mule_orientation = theta_list[0]
        self.end_trailer_pose = offset

        self.control_actions = []

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
        If $(x_{trailer, N}, y_{trailer, N})$ is position of last trailer of the robot, $k^{th}$ hitch and trailer position is calculated as:
            (assuming hitch-trailers are numbered $1$ to $N$, from mule to last trailer, and there are total $N$ hitch-trailer pairs)

            $x_{trailer, k} = \sum_{i=k}^{N} D_{i}\cos(\phi_{i}) + L_{i}\cos(\theta_{i})$ and $y_{trailer, k} = \sum_{i=k}^{N} D_{i}\sin(\phi_{i}) + L_{i}\sin(\theta_{i})$

            $x_{hitch, k} = x_{trailer, k} + D_{k-1}\cos(\phi_{k-1})$ and $y_{hitch, k} = y_{trailer, k} + D_{k-1}\sin(\phi_{k-1})$
        '''
        xi_trailer, yi_trailer = self.end_trailer_pose
        xi_hitch, yi_hitch = 0, 0

        try:
            for i in range(self.trailer_count - 1, k-1, -1):
                xi_trailer += self.D_list[i] * np.cos(self.phi_list[i]) + self.L_list[i] * np.cos(self.theta_list[i])
                yi_trailer += self.D_list[i] * np.sin(self.phi_list[i]) + self.L_list[i] * np.sin(self.theta_list[i])

            xi_hitch = xi_trailer + self.L_list[k-1] * np.cos(self.theta_list[k-1])
            yi_hitch = yi_trailer + self.L_list[k-1] * np.sin(self.theta_list[k-1])

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

        v_ht , omega_ht = v, omega
        for i in range(self.trailer_count):
            _, trailer_old = self.calculate_kth_hitch_trailer_pose(i+1)
            x_trailer_new = trailer_old[0] + self.v_trailer(v_ht, omega_ht, i) * np.cos(self.phi_list[i]) * dt
            y_trailer_new = trailer_old[1] + self.v_trailer(v_ht, omega_ht, i) * np.sin(self.phi_list[i]) * dt

            if i == self.trailer_count - 1:
                self.end_trailer_pose = [x_trailer_new, y_trailer_new]

            self.control_actions.append([v_ht, omega_ht])
            v_ht = self.v_trailer(v_ht, omega_ht, i)

            if i > 1:
                omega_ht = self.omega_next(v_ht, self.phi_dot(v_ht, omega_ht, i), i)

            self.theta_list[i] += omega_ht * dt
            self.phi_list[i] += self.phi_dot(v_ht, omega_ht, i) * dt

        self.mule_position = [x_mule_new, y_mule_new]
        self.mule_orientation = self.theta_list[0]

        return self.mule_position

        # try:
        #     print(self.mule_position, self.calculate_mule_pose())
        #     if all(self.mule_position[i] == self.calculate_mule_pose()[i] for i in range(self.trailer_count)):
        #         return self.mule_position, self.mule_orientation
        #     else:
        #         print("MISMATCH")
        #     print("MISMATCH")
        # except:
        #     print("Internal Error")
        #     return None

    def v_trailer(self, v, omega, i):
        r'''
        $v_{t, i} = v\cos(\theta_{i}-\phi_{i}) + L_{i}\omega\sin(\theta_{i}-\phi_{i})$
        '''
        return v * np.cos(self.phi_list[i] - self.theta_list[i]) + self.L_list[i] * omega * np.sin(self.phi_list[i] - self.theta_list[i])

    def phi_dot(self, v, omega, i):
        r'''
        $\dot\phi_{t, i} = \frac{v\sin(\theta_{i}-\phi_{i}) - L_{i}\omega\cos(\theta_{i}-\phi_{i})}{D_{i}}$
        '''
        return (v * np.sin(self.phi_list[i] - self.theta_list[i]) - self.L_list[i] * omega * np.cos(self.phi_list[i] - self.theta_list[i])) / self.D_list[i]

    def omega_next(self, v_ht, phi_dot, i):
        r'''
        $\dot\omega_{i} = \frac{v_{i} \sin(\phi_{i-1} - \theta_{i}) - D_{i-1} \dot\phi_{i-1} \cos(\phi_{i-1} - \theta_{i})}{L{i}}$
        '''
        return (v_ht * np.sin(self.phi_list[i] - self.theta_list[i+1]) - self.D_list[i] * phi_dot * np.cos(self.phi_list[i] - self.theta_list[i+1])) / self.L_list[i+1]

    def diagnostics(self):
        r'''
        Displays robot pose in tabular format
        +--------------+-------------+--------------+--------------+--------------+----....
        | Mule         | Hitch 1     | Trailer 1    | Hitch 2      | Trailer 2    |
        +==============+=============+==============+==============+==============+====....
        |[$x_{mule}$, $y_{mule}$] | [$x_{h, 1}$, $y_{h,1}$]  | [$x_{t, 1}$, $y_{t,1}$]   | [$x_{h, 2}$, $y_{h,2}$]   | [$x_{t, 2}$, $y_{t,2}$]   |
        +--------------+-------------+--------------+--------------+--------------+----....
        '''
        x_err = np.square(np.abs((self.mule_position[0] - self.calculate_mule_pose()[0]) / self.mule_position[0]))
        y_err = np.square(np.abs((self.mule_position[1] - self.calculate_mule_pose()[1]) / self.mule_position[1]))
        error = 100 * np.sqrt(x_err + y_err)

        headers = ["Mule"]
        rows = [[[round(float(xy), 2) for xy in self.calculate_mule_pose()]]]
        for i in range(self.trailer_count):
            headers += [f"Hitch {i+1}"] + [f"Trailer {i+1}"]
            hitch, trailer = self.calculate_kth_hitch_trailer_pose(i+1)
            hitch = [round(float(xy), 2) for xy in hitch]
            trailer = [round(float(xy), 2) for xy in trailer]
            rows[0].append(hitch)
            rows[0].append(trailer)
        print(tabulate(rows, headers, tablefmt='grid'))

        if error > 1 and error != 100:
            print(f"Percentage Error in Mule pose: {error:.2f}%. Precautionary measures advised.")


rd = RobotDynamics([0,0], [1.5, 1.5], [2.0, 2.0], [np.pi/4, np.pi/4], [np.pi/4, np.pi/4], trailer_count=2, direction=True)
rd.diagnostics()
for i in range(2):
    rd.update_state([1, 0.01], 0.1)
rd.diagnostics()

class Path():
    '''
    Generates parametric paths. Also returns the path's curvature.
    '''
    @staticmethod
    def circle(t, R=10.0, v=1.0):
        r'''
        Generates a circle path and its kinematic properties.
        Returns: $x$, $y$, $\theta$ (heading), $\kappa$ (curvature)
        '''
        omega = v / R
        x = R * np.cos(omega * t)
        y = R * np.sin(omega * t)
        theta = np.arctan2(y, x) + np.pi / 2
        kappa = 1.0 / R
        return x, y, theta, kappa

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
        return x, y, theta, kappa

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
        return x, y, theta, kappa

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
        return x, y, theta, kappa


class FeedbackController:
    pass

class MPC(FeedbackController):
    pass

class PID(FeedbackController):
    pass

class PP(FeedbackController):
    pass

class SMC(FeedbackController):
    pass

class RobotDynamics():
    pass


class Obstacle():
    pass


class PathPlanningFrenetFrame:
    def __init__(self, T, dt):
        self.T = T
        self.dt = dt

        self.current_states = []
        self.target_states = []
        self.control_actions = []

        self.current_states_frenet = []
        self.target_states_frenet = []
        self.control_actions_frenet = []

        self.control_loop()
        self.plot()

    def control_loop(self):
        robot = RobotDynamics([0,0], [1.5], [2.0], [np.pi/4], [np.pi/4], trailer_count=1, direction=True)
        current_state = robot.mule_position
        for t in np.arange(0, self.T, self.dt):
            target_state = Path.circle(t)
            control_action = FeedbackController.mpc(current_state, target_state)
            updated_state = robot.update_state(control_action)
            current_state = updated_state
            self.store(current_state, target_state, control_action)

    def store(self, current_state, target_state, control_action):
        self.current_states.append(current_state)
        self.target_states.append(target_state)
        self.control_actions.append(control_action)

    def plot(self):
        pass


# trajectory = PathPlanningFrenetFrame(60, 0.1)
