/*
    Header for path planning and task solver class
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#ifndef PATH_PLANNING_FRENET_FRAME_HPP
#define PATH_PLANNING_FRENET_FRAME_HPP

#include <cmath>
#include <array>
#include <chrono>
#include <vector>
#include <ostream>
// #include <matplotlibcpp.h>

#include "Path.hpp"
#include "Controller.hpp"
#include "RobotDynamics.hpp"

using namespace std;
// using namespace plt = matplotlibcpp;

class PathPlanningFrenetFrame {
    private:
        RobotDynamics& robot;
        Path& target_path;
        Controller& controller;
        float T;
        float dt;

        vector<PathPoint> current_states;
        vector<PathPoint> target_states;
        vector<vector<float>> control_actions;

        vector<vector<vector<float>>> hitches;
        vector<vector<vector<float>>> trailers;

        vector<vector<float>> current_states_frenet;
        vector<vector<float>> target_states_frenet;
        vector<vector<float>> control_actions_frenet;

    public:
        PathPlanningFrenetFrame(RobotDynamics&, Path&, Controller&, float, float);
        ~PathPlanningFrenetFrame();
        void store_hitch_trailer();
        void control_loop();
        void diagnostics();
        void plot(int);
        void display_time(float, float);
};

#endif
