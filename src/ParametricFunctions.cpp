/*
    Script for conversions between Cartesian and Frenet-Serret frame
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#include "ParametricFunctions.hpp"

ParametricFunctions::ParametricFunctions(double radius){
    this->radius = radius;
}

ParametricFunctions::ParametricFunctions(double major_axis, double minor_axis){
    this->minor_axis = minor_axis;
    this->major_axis = major_axis;
}

ParametricFunctions::~ParametricFunctions(){}

array<double,2> ParametricFunctions::figure_eight(double t) {
    double x = this->radius * sin(t);
    double y = this->radius * sin(t) * cos(t);
    return {x, y};
}

array<double,2> ParametricFunctions::figure_circle(double t) {
    double x = this->radius * sin(t);
    double y = this->radius * sin(t) * cos(t);
    return {x, y};
}

array<double,2> ParametricFunctions::figure_ellipse(double t) {
    double x = this->major_axis * cos(t);
    double y = this->minor_axis * sin(t);
    return {x, y};
}
