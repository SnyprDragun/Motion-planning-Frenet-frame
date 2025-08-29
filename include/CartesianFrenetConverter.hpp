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
#include <stdexcept>

using namespace std;

class CartesianFrenetConverter {
    private:
        

    public:
        CartesianFrenetConverter();
        ~CartesianFrenetConverter();
        static pair<array<double, 3>, array<double, 3>> cartesian_to_frenet(double, double, double, double, double, double, double, double, double, double, double, double);
        static array<double, 6> frenet_to_cartesian(double, double, double, double, double, double, array<double, 3>, array<double, 3>);
        static double normalize_angle(double);
};

#endif
