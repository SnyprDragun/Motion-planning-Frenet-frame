/*
    Header for conversions between Cartesian and Frenet-Serret frame
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#ifndef CARTESIAN_FRENET_CONVERTER_HPP
#define CARTESIAN_FRENET_CONVERTER_HPP

#include <array>
#include <cmath>
#include <utility>
#include <iostream>
#include <stdexcept>

using namespace std;

class CartesianFrenetConverter {
    private:
        

    public:
        CartesianFrenetConverter();
        ~CartesianFrenetConverter();

        static pair<array<double, 3>, array<double, 3>> cartesian_to_frenet(
            double rs, 
            double rx, 
            double ry, 
            double rtheta,
            double rkappa, 
            double rdkappa,
            double x, 
            double y, 
            double v, 
            double a,
            double theta, 
            double kappa
        );

        static array<double, 6> frenet_to_cartesian(
            double rs, 
            double rx, 
            double ry, 
            double rtheta,
            double rkappa, 
            double rdkappa,
            array<double, 3> s_condition,
            array<double, 3> d_condition
        );

        static double normalize_angle(double angle);
};

#endif
