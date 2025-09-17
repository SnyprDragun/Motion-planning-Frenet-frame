/*
    Header to simulate dynamics of chosen robot
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#ifndef ROBOT_DYNAMICS_HPP
#define ROBOT_DYNAMICS_HPP

#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <string>

using namespace std;

class RobotDynamics {
    private:
        

    public:
        array<double,2> end_trailer_pose;
        vector<double> L_list;
        vector<double> D_list;
        vector<double> theta_list;
        vector<double> phi_list;
        array<double,2> mule_position;
        vector<array<double,2>> control_actions;
        int trailer_count;
        double mule_orientation;
        bool direction;

        RobotDynamics(array<double,2>, vector<double>, vector<double>, vector<double>, vector<double>, int, bool);
        array<double,2> calculate_mule_pose();
        pair<array<double,2>, array<double,2>> calculate_kth_hitch_trailer_pose(int);
        array<double,3> update_state(const array<double,2>, double);
        double v_trailer(double, double, int);
        double phi_dot(double, double, int);
        double omega_next(double, double, int);
        void diagnostics();
};

#endif
