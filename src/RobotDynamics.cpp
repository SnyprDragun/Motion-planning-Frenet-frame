/*
    Script to simulate dynamics of chosen robot
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#include "RobotDynamics.hpp"

/*
Dynamics are given as follows:
    > The position of the robot (mule) is given by $(x, y, \theta)$. This point is at the centre of the front axle. 
    > Along the line from this point to the rear wheel and beyond is the hitch, at distance $L$.
    > Further, we have a trailer connected to this hitch point. 
    > The rear wheel of this trailer is at a distance $D$ from the hitch. 
    > We assume the line connecting the hitch and the two back wheels form an angle $\phi$.
*/

RobotDynamics::RobotDynamics(array<double,2> offset, vector<double> L_list, vector<double> D_list, vector<double> theta_list, vector<double> phi_list, int trailer_count, bool direction){
    this->end_trailer_pose = offset;
    this->L_list = L_list;
    this->D_list = D_list;
    this->theta_list = theta_list;
    this->phi_list = phi_list;
    this->trailer_count = trailer_count;
    this->direction = direction;

    this->mule_position = {0, 0};
    this->mule_orientation = this->theta_list.empty() ? 0.0 : this->theta_list[0];
    this->control_actions = {};

    calculate_mule_pose();
}

array<double,2> RobotDynamics::calculate_mule_pose(){
    /*
    If $(x_{trailer, N}, y_{trailer, N})$ is position of last trailer of the robot (offset), then mule position is calculated as:

        $x_{mule} = x_{trailer, N} + \sum_{i=1}^{N} D_{i}\cos(\phi_{i}) + L_{i}\cos(\theta_{i})$ and $y_{mule} = y_{trailer,N} + \sum_{i=1}^{N} D_{i}\sin(\phi_{i}) + L_{i}\sin(\theta_{i})$ where $N$ represents total numver of trailers
    */

    array<double,2> mule_position = {0.0, 0.0};

    for (int i = 0; i < trailer_count; i++) {
        double xi = this->D_list[i] * cos(this->phi_list[i]) + this->L_list[i] * cos(this->theta_list[i]);
        double yi = this->D_list[i] * sin(this->phi_list[i]) + this->L_list[i] * sin(this->theta_list[i]);

        mule_position[0] += xi;
        mule_position[1] += yi;
    }

    mule_position[0] += this->end_trailer_pose[0];
    mule_position[1] += this->end_trailer_pose[1];
    this->mule_position = mule_position;

    return this->mule_position;
}

pair<array<double,2>, array<double,2>> RobotDynamics::calculate_kth_hitch_trailer_pose(int k){
    /*
    If $(x_{mule}, y_{mule})$ is position of the mule, $k^{th}$ hitch and trailer position is calculated as:
        (assuming hitch-trailers are numbered $1$ to $N$, from mule to last trailer, and there are total $N$ hitch-trailer pairs)

        $x_{trailer, k} = x_{mule} - \sum_{i=1}^{k} D_{i}\cos(\phi_{i}) + L_{i}\cos(\theta_{i})$ and $y_{trailer, k} = y_{mule} - \sum_{i=1}^{k} D_{i}\sin(\phi_{i}) + L_{i}\sin(\theta_{i})$

        $x_{hitch, k} = x_{trailer, k} + D_{k}\cos(\phi_{k})$ and $y_{hitch, k} = y_{trailer, k} + D_{k}\sin(\phi_{k})$
    */

    double xi_trailer = this->mule_position[0];
    double yi_trailer = this->mule_position[1];
    double xi_hitch = 0.0, yi_hitch = 0.0;

    if (k <= 0 || k > this->trailer_count) {
        cerr << "ERROR: Robot has only " << this->trailer_count 
             << " hitch-trailers, you are asking pose for " << k << endl;
        cerr << "Proceeding with last hitch-trailer pose..." << endl;
        return calculate_kth_hitch_trailer_pose(this->trailer_count);
    }

    for (int i = 0; i < k; i++) {
        xi_trailer -= this->D_list[i] * cos(this->phi_list[i]) + this->L_list[i] * cos(this->theta_list[i]);
        yi_trailer -= this->D_list[i] * sin(this->phi_list[i]) + this->L_list[i] * sin(this->theta_list[i]);
    }

    xi_hitch = xi_trailer + this->D_list[k-1] * cos(this->phi_list[k-1]);
    yi_hitch = yi_trailer + this->D_list[k-1] * sin(this->phi_list[k-1]);

    array<double,2> hitch = {xi_hitch, yi_hitch};
    array<double,2> trailer = {xi_trailer, yi_trailer};

    return {hitch, trailer};
}

array<double,3> RobotDynamics::update_state(const array<double,2> control_action, double dt){
    /*
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
    */

    double v = control_action[0];
    double omega = control_action[1];

    double x_mule_new = mule_position[0] + v * cos(theta_list[0]) * dt;
    double y_mule_new = mule_position[1] + v * sin(theta_list[0]) * dt;

    vector<double> new_theta_list = theta_list;
    vector<double> new_phi_list = phi_list;

    double v_ht = v;
    double omega_ht = omega;
    double v_ht_new = 0.0, omega_ht_new = 0.0;

    for (int i = 0; i < trailer_count; i++) {
        control_actions.push_back({v_ht, omega_ht});

        auto trailer_pair = calculate_kth_hitch_trailer_pose(i+1);
        array<double,2> trailer_old = trailer_pair.second;

        double x_trailer_new = trailer_old[0] + v_trailer(v_ht, omega_ht, i) * cos(phi_list[i]) * dt;
        double y_trailer_new = trailer_old[1] + v_trailer(v_ht, omega_ht, i) * sin(phi_list[i]) * dt;

        if (i == trailer_count - 1) {
            end_trailer_pose = {x_trailer_new, y_trailer_new};
        }

        if (i > 0) {
            omega_ht_new = omega_next(v_ht, phi_dot(v_ht, omega_ht, i), i);
            v_ht_new = v_trailer(v_ht, omega_ht, i);
            omega_ht = omega_ht_new;
            v_ht = v_ht_new;
        }

        new_theta_list[i] += omega_ht * dt;
        new_phi_list[i] += phi_dot(v_ht, omega_ht, i) * dt;
    }

    theta_list = new_theta_list;
    phi_list = new_phi_list;
    mule_position = {x_mule_new, y_mule_new};
    mule_orientation = theta_list[0];

    return {mule_position[0], mule_position[1], mule_orientation};
}

double RobotDynamics::v_trailer(double v, double omega, int i){
    /*
    $v_{t, i} = v\cos(\theta_{i}-\phi_{i}) + L_{i}\omega\sin(\theta_{i}-\phi_{i})$
    */

    return v * cos(theta_list[i] - phi_list[i]) + L_list[i] * omega * sin(theta_list[i] - phi_list[i]);
}

double RobotDynamics::phi_dot(double v, double omega, int i){
    /*
    $\dot\phi_{t, i} = \frac{v\sin(\theta_{i}-\phi_{i}) - L_{i}\omega\cos(\theta_{i}-\phi_{i})}{D_{i}}$
    */

    return (v * sin(theta_list[i] - phi_list[i]) - L_list[i] * omega * cos(theta_list[i] - phi_list[i])) / D_list[i];
}

double RobotDynamics::omega_next(double v_ht, double phi_dot, int i){
    /*
    $\omega_{i} = \frac{v_{hitch, i-1} \sin(\phi_{i-1} - \theta_{i}) - D_{i-1} \dot\phi_{i-1} \cos(\phi_{i-1} - \theta_{i})}{L{i}}$
    */

    return (v_ht * sin(phi_list[i-1] - theta_list[i]) - D_list[i-1] * phi_dot * cos(phi_list[i-1] - theta_list[i])) / L_list[i];
}

void RobotDynamics::diagnostics(){
    /*
    Displays robot pose in tabular format
    +--------------+-------------+--------------+--------------+--------------+----....
    | Mule         | Hitch 1     | Trailer 1    | Hitch 2      | Trailer 2    |
    +==============+=============+==============+==============+==============+====....
    |[$x_{mule}$, $y_{mule}$] | [$x_{h, 1}$, $y_{h,1}$]  | [$x_{t, 1}$, $y_{t,1}$]   | [$x_{h, 2}$, $y_{h,2}$]   | [$x_{t, 2}$, $y_{t,2}$]   |
    +--------------+-------------+--------------+--------------+--------------+----....
    */

    cout << setw(12) << "Mule";
    for (int i = 0; i < trailer_count; i++) {
        cout << setw(12) << ("Hitch " + to_string(i+1));
        cout << setw(12) << ("Trailer " + to_string(i+1));
    }
    cout << "\n";

    cout << fixed << setprecision(2);
    cout << setw(12) << ("(" + to_string(mule_position[0]) + "," + to_string(mule_position[1]) + ")");

    for (int i = 0; i < trailer_count; i++) {
        auto poses = calculate_kth_hitch_trailer_pose(i+1);
        auto hitch = poses.first;
        auto trailer = poses.second;

        cout << setw(12) << ("(" + to_string(hitch[0]) + "," + to_string(hitch[1]) + ")");
        cout << setw(12) << ("(" + to_string(trailer[0]) + "," + to_string(trailer[1]) + ")");
    }
    cout << "\n";
}
