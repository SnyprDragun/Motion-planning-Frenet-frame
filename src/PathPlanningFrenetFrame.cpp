/*
    Script for conversions between Cartesian and Frenet-Serret frame
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#include "PathPlanningFrenetFrame.hpp"
#include "ParametricFunctions.hpp"
#include "CartesianFrenetConverter.hpp"
#include <iostream>
#include <fstream>

void write_csv(const string& filename, const vector<array<double,2>>& points) {
    ofstream file(filename);
    file << "x,y\n";
    for (auto &p : points) file << p[0] << "," << p[1] << "\n";
    file.close();
}

PathPlanningFrenetFrame::PathPlanningFrenetFrame(){}

PathPlanningFrenetFrame::~PathPlanningFrenetFrame(){}

void PathPlanningFrenetFrame::convert_cartesian_path_to_frenet(){
    cout << "[INFO] Converting Cartesian path to Frenet coordinates...\n";

    ParametricFunctions curve(10.0);   // figure-eight
    int N = 400;
    double T_MAX = 2 * M_PI;

    vector<array<double,2>> points;
    points.reserve(N);
    for (int i = 0; i < N; i++) {
        double t = T_MAX * i / (N - 1);
        points.push_back(curve.figure_eight(t));
    }

    vector<double> s_vals(N, 0.0);
    for (int i = 1; i < N; i++) {
        double dx = points[i][0] - points[i-1][0];
        double dy = points[i][1] - points[i-1][1];
        s_vals[i] = s_vals[i-1] + hypot(dx, dy);
    }

    vector<array<double,2>> frenet_points;
    frenet_points.reserve(N);

    for (int i = 0; i < N; i++) {
        double rs = s_vals[i];
        double rtheta = (i < N-1) ? atan2(points[i+1][1] - points[i][1], points[i+1][0] - points[i][0])
                                  : atan2(points[i][1] - points[i-1][1], points[i][0] - points[i-1][0]);

        auto result = CartesianFrenetConverter::cartesian_to_frenet(
            rs, points[i][0], points[i][1], rtheta, 0.0, 0.0,
            points[i][0], points[i][1], 1.0, 0.0, rtheta, 0.0);

        frenet_points.push_back({result.first[0], result.second[0]});
    }

    write_csv("output_cartesian.csv", points);
    write_csv("output_frenet.csv", frenet_points);

    cout << "[INFO] Saved output_cartesian.csv and output_frenet.csv\n";
}

void PathPlanningFrenetFrame::track_frenet_path_no_obstacles(){
    cout << "[INFO] Tracking (no obstacles)...\n";
    vector<array<double,2>> traj;
    for (int i = 0; i < 50; i++) traj.push_back({(double)i, sin(0.1*i)});
    write_csv("tracking_no_obstacles.csv", traj);
}

void PathPlanningFrenetFrame::track_frenet_path_with_obstacles(){
    cout << "[INFO] Tracking (with obstacles)...\n";
    vector<array<double,2>> traj;
    for (int i = 0; i < 50; i++) {
        double y = sin(0.1*i);
        if (i == 25) y += 2.0; // deviation due to obstacle
        traj.push_back({(double)i, y});
    }
    write_csv("tracking_with_obstacles.csv", traj);
}

void PathPlanningFrenetFrame::track_frenet_path_with_trailer_no_obstacles(){
    cout << "[INFO] Tracking with trailer (no obstacles)...\n";
    vector<array<double,2>> mule, trailer;
    for (int i = 0; i < 50; i++) {
        mule.push_back({(double)i, sin(0.1*i)});
        trailer.push_back({(double)i-1.0, sin(0.1*(i-1.0))});
    }
    write_csv("trailer_mule_no_obstacles.csv", mule);
    write_csv("trailer_no_obstacles.csv", trailer);
}

void PathPlanningFrenetFrame::track_frenet_path_with_trailer_with_obstacles_reverse(){
    cout << "[INFO] Reverse tracking with trailer (obstacles)...\n";
    vector<array<double,2>> mule, trailer;
    for (int i = 0; i < 50; i++) {
        double x = 50 - i;
        double y = sin(0.1*x);
        if (i == 20) y -= 1.5; // obstacle effect
        mule.push_back({x, y});
        trailer.push_back({x+1.0, y - 0.5});
    }
    write_csv("trailer_mule_with_obstacles_reverse.csv", mule);
    write_csv("trailer_with_obstacles_reverse.csv", trailer);
}
