/*
    Main script
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#include <chrono>
#include "RobotDynamics.hpp"

using namespace chrono;

int main() {
    auto start = high_resolution_clock::now();

    array<double,2> OFFSET = {10.0, -7.0};
    vector<double> THETA = {M_PI/2, M_PI/2};
    vector<double> PHI   = {M_PI/2, M_PI/2};
    vector<double> L     = {1.5, 1.5};
    vector<double> D     = {2.0, 2.0};

    RobotDynamics robot(OFFSET, L, D, THETA, PHI, 2, true);
    robot.diagnostics();

    auto end = high_resolution_clock::now();
    duration<double> elapsed = end - start;
    cout << "Execution time: " << elapsed.count() << " seconds\n";

    return 0;
}
