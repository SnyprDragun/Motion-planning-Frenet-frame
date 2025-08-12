/*
    Script for conversions between Cartesian and Frenet-Serret frame
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#include "ParametricFunctions.hpp"

ParametricFunctions::ParametricFunctions(){}

ParametricFunctions::~ParametricFunctions(){}

PathPoint ParametricFunctions::straight_line(double t, double slope=0.0, double intercept=0.0, double v=1.0) {
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

PathPoint ParametricFunctions::circle(double t, double R=10.0, double v=1.0) {
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

PathPoint ParametricFunctions::ellipse(double t, double a=10.0, double b=6.0, double v=1.0) {
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

PathPoint ParametricFunctions::figure_eight(double t, double a=10.0, double v=1.0) {
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
