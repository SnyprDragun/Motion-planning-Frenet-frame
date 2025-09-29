#!/Users/subhodeep/venv/bin/python3
'''class to simulate dynamics of chosen robot'''
import numpy as np
from tabulate import tabulate

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
