/*
    Main script
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame

    To compile:
        g++ -std=c++17 \
        -Iinclude \
        -I/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/ \
        -I/opt/homebrew/include/ \
        -I/opt/homebrew/opt/python@3.13/Frameworks/Python.framework/Versions/3.13/include/python3.13 \
        -I/Users/subhodeep/venv/lib/python3.13/site-packages/numpy/_core/include \
        src/main.cpp src/RobotDynamics.cpp src/Controller.cpp src/PathPlanningFrenetFrame.cpp src/Path.cpp \
        -L/opt/homebrew/lib \
        -L/opt/homebrew/opt/python@3.13/Frameworks/Python.framework/Versions/3.13/lib \
        -lnlopt \
        -lpython3.13 \
        -Wno-deprecated-declarations \
        -o main

    To run:
        ./main
*/

#include "PathPlanningFrenetFrame.hpp"

int main() {
    auto start = chrono::high_resolution_clock::now();

    vector<float> OFFSET = {10.0, -7.0};
    vector<float> THETA = {M_PI/2, M_PI/2};
    vector<float> PHI   = {M_PI/2, M_PI/2};
    vector<float> L     = {1.5, 1.5};
    vector<float> D     = {2.0, 2.0};

    RobotDynamics robot(OFFSET, L, D, THETA, PHI, 2, true);
    robot.diagnostics();
    Path target_path("circle", 10.0, 12.0, 8.0, 0.0, 0.0, 1.0);
    MPCParams params(0.1, 20);   // dt=0.1, N=20
    Controller mpc(params);

    PathPlanningFrenetFrame trajectory(robot, target_path, mpc, 10.0, 0.1);
    trajectory.control_loop();
    trajectory.diagnostics();
    robot.diagnostics();

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<float> elapsed = end - start;
    cout << "Execution time: " << elapsed.count() << " seconds\n";

    trajectory.plot(100);

    return 0;
}
