/*
    Header for point coordinate on parametric path
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#ifndef PATH_POINT_HPP
#define PATH_POINT_HPP

struct PathPoint {
    double x;
    double y;
    double theta;
    double kappa;
};

#endif
