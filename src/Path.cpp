/*
    Script for parametric equations of candidate path
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#include "Path.hpp"

Path::Path(string shape, float R=10.0, float a=12.0, float b=8.0, float slope=0.0, float intercept=0.0, float v=1.0){
    this->shape = shape;
    this->R = R;
    this->a = a;
    this->b = b;
    this->slope = slope;
    this->intercept = intercept;
    this->v = v;
    this->obstacles = 0;
}

Path::~Path(){}

PathPoint Path::straight_line(double t, double slope=0.0, double intercept=0.0, double v=1.0) {
    /*
    Generates a straight line path and its kinematic properties.
    Returns: x, y, theta (heading), kappa (curvature)
    */

    PathPoint path_point;

    path_point.x = v * t;
    path_point.y = slope * path_point.x + intercept;
    path_point.theta = atan2(slope, 1.0);
    path_point.kappa = 0.0;

    return path_point;
}

PathPoint Path::circle(double t, double R=10.0, double v=1.0) {
    /*
    Generates a circle path and its kinematic properties.
    Returns: x, y, theta (heading), kappa (curvature)
    */

    PathPoint path_point;

    double omega = v / R;
    path_point.x = R * cos(omega * t);
    path_point.y = R * sin(omega * t);
    path_point.theta = atan2(path_point.y, path_point.x) + M_PI / 2;
    path_point.kappa = 1.0 / R;

    return path_point;
}

PathPoint Path::ellipse(double t, double a=10.0, double b=6.0, double v=1.0) {
    /*
    Generates an ellipse path and its kinematic properties.
    Returns: x, y, theta (heading), kappa (curvature)
    */

    PathPoint path_point;

    double omega = v / ((a + b) / 2.0);

    path_point.x = a * cos(omega * t);
    path_point.y = a * sin(omega * t);

    double dx = -a * omega * sin(omega * t);
    double dy = b * omega * cos(omega * t);
    double ddx = -a * omega * omega * cos(omega * t);
    double ddy = -b * omega * omega * sin(omega * t);

    path_point.theta = atan2(dy, dx);

    double speed = dx * dx + dy * dy;
    if (speed < 1e-6) {
        path_point.kappa = 0;
    }
    else{
        path_point.kappa = (dx * ddy - dy * ddx) / pow(speed, 1.5);
    }

    return path_point;
}

PathPoint Path::figure_eight(double t, double a=10.0, double v=1.0) {
    /*
    Generates a figure-eight path and its kinematic properties.
    Returns: x, y, theta (heading), kappa (curvature)
    */

    PathPoint path_point;

    double omega = v / a;

    double dx = a * omega * cos(omega * t);
    double dy = a * omega * sin(2 * omega * t);
    double ddx = -a * omega * omega * sin(omega * t);
    double ddy = -2 * a * omega * omega * sin(2 * omega * t);

    path_point.x = a * sin(omega * t);
    path_point.y = a * sin(omega * t) * cos(omega * t);
    path_point.theta = atan2(dy, dx);

    double speed = dx * dx + dy * dy;
    if (speed < 1e-6) {
        path_point.kappa = 0;
    }
    else{
        path_point.kappa = (dx * ddy - dy * ddx) / pow(speed, 1.5);
    }

    return path_point;
}

PathPoint Path::equation(float t){
    try {
        if (this->shape == "circle") {
            return Path::circle(t, R, v);
        } else if (shape == "ellipse") {
            return Path::ellipse(t, a, b, v);
        } else if (shape == "figure_eight") {
            return Path::figure_eight(t, a, v);
        } else if (shape == "straight_line") {
            return Path::straight_line(t, slope, intercept, v);
        } else {
            throw invalid_argument("Unknown path shape!");
        }
    } catch (const exception& e) {
        cerr << "Provide correct path format! (" << e.what() << ")\n";
        return {0.0, 0.0}; // fallback
    }
}

void Path::add_obstacle(){}
