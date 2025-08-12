/*
    Header for conversions between Cartesian and Frenet-Serret frame
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#ifndef PARAMETRIC_FUNCTIONS_HPP
#define PARAMETRIC_FUNCTIONS_HPP

#include <array>
#include <cmath>

using namespace std;

struct PathPoint {
    double x;
    double y;
    double theta;
    double kappa;
};

class ParametricFunctions {
    private:
        

    public:
        ParametricFunctions();
        ~ParametricFunctions();
        static PathPoint straight_line(double, double, double, double);
        static PathPoint circle(double, double, double);
        static PathPoint ellipse(double, double, double, double);
        static PathPoint figure_eight(double, double, double);
};

#endif
