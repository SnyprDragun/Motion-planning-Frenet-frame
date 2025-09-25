/*
    Main script
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#include <chrono>
#include "PathPlanningFrenetFrame.hpp"

using namespace chrono;

int main() {
    auto start = high_resolution_clock::now();

    vector<float> OFFSET = {10.0, -7.0};
    vector<float> THETA = {M_PI/2, M_PI/2};
    vector<float> PHI   = {M_PI/2, M_PI/2};
    vector<float> L     = {1.5, 1.5};
    vector<float> D     = {2.0, 2.0};

    RobotDynamics robot(OFFSET, L, D, THETA, PHI, 2, true);
    robot.diagnostics();
    Path target_path("circle", 10.0, 12.0, 8.0, 0.0, 0.0, 1.0);
    Controller controller(1.0, 0.1, 0.01);

    PathPlanningFrenetFrame path_planner(robot, target_path, controller, 10.0, 0.1);

    auto end = high_resolution_clock::now();
    duration<float> elapsed = end - start;
    cout << "Execution time: " << elapsed.count() << " seconds\n";

    return 0;
}
