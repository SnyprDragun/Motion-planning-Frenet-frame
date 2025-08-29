/*
    Script for conversions between Cartesian and Frenet-Serret frame
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#include "CartesianFrenetConverter.hpp"

CartesianFrenetConverter::CartesianFrenetConverter(){}

CartesianFrenetConverter::~CartesianFrenetConverter(){}

pair<array<double, 3>, array<double, 3>> CartesianFrenetConverter::cartesian_to_frenet(double rs, double rx, double ry, double rtheta, double rkappa, double rdkappa, double x, double y, double v, double a, double theta, double kappa){
    /*
    Convert state from Cartesian coordinate to Frenet coordinate

    Parameters
    ----------
        rs: reference line s-coordinate
        rx, ry: reference point coordinates
        rtheta: reference point heading
        rkappa: reference point curvature
        rdkappa: reference point curvature rate
        x, y: current position
        v: velocity
        a: acceleration
        theta: heading angle
        kappa: curvature

    Returns
    -------
        s_condition: [s(t), s'(t), s''(t)]
        d_condition: [d(s), d'(s), d''(s)]
    */

    double dx = x - rx;
    double dy = y - ry;

    double cos_theta_r = cos(rtheta);
    double sin_theta_r = sin(rtheta);

    double cross_rd_nd = cos_theta_r * dy - sin_theta_r * dx;
    double d = copysign(hypot(dx, dy), cross_rd_nd);

    double delta_theta = theta - rtheta;
    double tan_delta_theta = tan(delta_theta);
    double cos_delta_theta = cos(delta_theta);

    double one_minus_kappa_r_d = 1 - rkappa * d;
    double d_dot = one_minus_kappa_r_d * tan_delta_theta;

    double kappa_r_d_prime = rdkappa * d + rkappa * d_dot;

    double d_ddot = (- kappa_r_d_prime * tan_delta_theta + one_minus_kappa_r_d / (cos_delta_theta * cos_delta_theta) * (kappa * one_minus_kappa_r_d / cos_delta_theta - rkappa));

    double s = rs;
    double s_dot = v * cos_delta_theta / one_minus_kappa_r_d;

    double delta_theta_prime = one_minus_kappa_r_d / cos_delta_theta * kappa - rkappa;
    double s_ddot = (a * cos_delta_theta - s_dot * s_dot * (d_dot * delta_theta_prime - kappa_r_d_prime)) / one_minus_kappa_r_d;

    return {
        array<double, 3>{s, s_dot, s_ddot}, 
        array<double, 3>{d, d_dot, d_ddot}
    };
}

array<double, 6> CartesianFrenetConverter::frenet_to_cartesian(double rs, double rx, double ry, double rtheta, double rkappa, double rdkappa, array<double, 3> s_condition, array<double, 3> d_condition){
    /*
    Convert state from Frenet coordinate to Cartesian coordinate

    Parameters
    ----------
        rs: reference line s-coordinate
        rx, ry: reference point coordinates
        rtheta: reference point heading
        rkappa: reference point curvature
        rdkappa: reference point curvature rate
        s_condition: [s(t), s'(t), s''(t)]
        d_condition: [d(s), d'(s), d''(s)]

    Returns
    -------
        x, y: position
        theta: heading angle
        kappa: curvature
        v: velocity
        a: acceleration
    */

    if (fabs(rs - s_condition[0]) >= 1.0e-6){
        throw std::runtime_error("The reference point s and s_condition[0] don't match");
    }

    double cos_theta_r = cos(rtheta);
    double sin_theta_r = sin(rtheta);

    double x = rx - sin_theta_r * d_condition[0];
    double y = ry + cos_theta_r * d_condition[0];

    double one_minus_kappa_r_d = 1 - rkappa * d_condition[0];

    double tan_delta_theta = d_condition[1] / one_minus_kappa_r_d;
    double delta_theta = atan2(d_condition[1], one_minus_kappa_r_d);
    double cos_delta_theta = cos(delta_theta);

    double theta = normalize_angle(delta_theta + rtheta);

    double kappa_r_d_prime = rdkappa * d_condition[0] + rkappa * d_condition[1];

    double kappa = (((d_condition[2] + kappa_r_d_prime * tan_delta_theta) * cos_delta_theta * cos_delta_theta) / one_minus_kappa_r_d + rkappa) * cos_delta_theta / one_minus_kappa_r_d;

    double d_dot = d_condition[1] * s_condition[1];
    double v = sqrt(one_minus_kappa_r_d * one_minus_kappa_r_d * s_condition[1] * s_condition[1] + d_dot * d_dot);

    double delta_theta_prime = one_minus_kappa_r_d / cos_delta_theta * kappa - rkappa;

    double a = (s_condition[2] * one_minus_kappa_r_d / cos_delta_theta + s_condition[1] * s_condition[1] / cos_delta_theta * (d_condition[1] * delta_theta_prime - kappa_r_d_prime));

    return array<double, 6>{x, y, theta, kappa, v, a};
}

double CartesianFrenetConverter::normalize_angle(double angle){
    /*
    Normalize angle to [-pi, pi]
    */

    const double two_pi = 2.0 * M_PI;
    return angle - two_pi * floor((angle + M_PI) / two_pi) - M_PI;
}
