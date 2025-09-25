/*
    Header for parametric equations of candidate path
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#ifndef PATH_HPP
#define PATH_HPP

#include <string>
#include <array>
#include <cmath>
#include <iostream>

#include "PathPoint.hpp"

using namespace std;

class Path {
    private:
        string shape;
        float R;
        float a;
        float b;
        float slope;
        float intercept;
        float v;
        float obstacles;

    public:
        Path(string, float, float, float, float, float, float);
        ~Path();
        static PathPoint straight_line(double, double, double, double);
        static PathPoint circle(double, double, double);
        static PathPoint ellipse(double, double, double, double);
        static PathPoint figure_eight(double, double, double);
        PathPoint equation(float);
        void add_obstacle();
};

#endif
