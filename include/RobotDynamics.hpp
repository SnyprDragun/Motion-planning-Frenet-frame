/*
    Header to simulate dynamics of chosen robot
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#ifndef ROBOT_DYNAMICS_HPP
#define ROBOT_DYNAMICS_HPP

#include <array>
#include <vector>
#include <string>
#include <iomanip>
#include <iostream>

#include "PathPoint.hpp"

using namespace std;

class RobotDynamics {
    private:
        

    public:
        vector<float> end_trailer_pose;
        vector<float> L_list;
        vector<float> D_list;
        vector<float> theta_list;
        vector<float> phi_list;
        vector<float> mule_position;
        vector<vector<float>> control_actions;
        int trailer_count;
        float mule_orientation;
        bool direction;

        RobotDynamics(vector<float>, vector<float>, vector<float>, vector<float>, vector<float>, int, bool);
        vector<float> calculate_mule_pose();
        pair<vector<float>, vector<float>> calculate_kth_hitch_trailer_pose(int);
        PathPoint update_state(const vector<float>, float);
        float v_trailer(float, float, int);
        float phi_dot(float, float, int);
        float omega_next(float, float, int);
        void diagnostics();
};

#endif
