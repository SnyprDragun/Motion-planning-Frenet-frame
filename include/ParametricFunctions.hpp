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

class ParametricFunctions {
    private:
        

    public:
        double radius;
        double minor_axis;
        double major_axis;
        ParametricFunctions(double);
        ParametricFunctions(double, double);
        ~ParametricFunctions();
        array<double,2> figure_eight(double);
        array<double,2> figure_circle(double);
        array<double,2> figure_ellipse(double);
};

#endif
