/*
    Script for controller class
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#include "Controller.hpp"

Controller::Controller(int N, double dt, double v,
                       double w_d_mule, double w_d_hitch, double w_d_trailer,
                       double w_d_dot_trailer, double w_phi,
                       double w_omega_mule, double w_omega_rate_mule,
                       pair<double, double> omega_bounds, double max_delta_omega,
                       bool use_exp_smooth, double smooth_alpha, bool is_reverse)
    : N(N), dt(dt), v(v),
      w_d_mule(w_d_mule), w_d_hitch(w_d_hitch), w_d_trailer(w_d_trailer),
      w_d_dot_trailer(w_d_dot_trailer), w_phi(w_phi),
      w_omega_mule(w_omega_mule), w_omega_rate_mule(w_omega_rate_mule),
      omega_min(omega_bounds.first), omega_max(omega_bounds.second),
      prev_omega(0.0), max_delta_omega(max_delta_omega),
      use_exp_smooth(use_exp_smooth), smooth_alpha(smooth_alpha),
      is_reverse(is_reverse)
{
    if (is_reverse) {
        this->w_d_trailer = 3500.0;
        this->w_omega_rate_mule = 9000.0;
    } else {
        this->w_d_trailer = 0.0;
        this->w_omega_rate_mule = 1.0;
    }
}

Controller::~Controller(){}

Eigen::Vector4d Controller::step_dynamics_forward(const Eigen::Vector4d& state, double omega, double L, double D) const {
    double x_mule = state[0];
    double y_mule = state[1];
    double theta_mule = state[2];
    double phi_hitch = state[3];

    double dx = v * cos(theta_mule);
    double dy = v * sin(theta_mule);
    double dtheta = omega;
    double dphi = (v * sin(theta_mule - phi_hitch) 
                  - L * omega * cos(theta_mule - phi_hitch)) / D;

    double x_mule_new = x_mule + dx * dt;
    double y_mule_new = y_mule + dy * dt;
    double theta_mule_new = CartesianFrenetConverter::normalize_angle(theta_mule + dtheta * dt);
    double phi_hitch_new = CartesianFrenetConverter::normalize_angle(phi_hitch + dphi * dt);

    return Eigen::Vector4d(x_mule_new, y_mule_new, theta_mule_new, phi_hitch_new);
}

Eigen::Vector4d Controller::step_dynamics_reverse(const Eigen::Vector4d& state, double omega, double L, double D) const {
    double x_trailer = state[0];
    double y_trailer = state[1];
    double theta_trailer = state[2];
    double phi_hitch = state[3];

    double theta_mule = CartesianFrenetConverter::normalize_angle(theta_trailer + phi_hitch);
    double vx_h = -v * cos(theta_mule) + L * omega * sin(theta_mule);
    double vy_h = -v * sin(theta_mule) - L * omega * cos(theta_mule);

    double vtrailer = vx_h * cos(theta_trailer) + vy_h * sin(theta_trailer);
    double vperp = -vx_h * sin(theta_trailer) + vy_h * cos(theta_trailer);

    double dtheta_trailer = vperp / D;
    double dphi = omega - dtheta_trailer;

    double x_trailer_new = x_trailer + vtrailer * cos(theta_trailer) * dt;
    double y_trailer_new = y_trailer + vtrailer * sin(theta_trailer) * dt;
    double theta_trailer_new = CartesianFrenetConverter::normalize_angle(theta_trailer + dtheta_trailer * dt);
    double phi_hitch_new = CartesianFrenetConverter::normalize_angle(phi_hitch + dphi * dt);

    return Eigen::Vector4d(x_trailer_new, y_trailer_new, theta_trailer_new, phi_hitch_new);
}

