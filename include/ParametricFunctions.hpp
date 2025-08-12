/*
    Header for parametric equations of candidate path
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#ifndef PARAMETRIC_FUNCTIONS_HPP
#define PARAMETRIC_FUNCTIONS_HPP

#include <array>
#include <cmath>
#include "PathPoint.hpp"

using namespace std;

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
