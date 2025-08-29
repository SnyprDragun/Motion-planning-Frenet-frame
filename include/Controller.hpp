/*
    Header for controller class
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#ifndef CONTROLLER_HPP
#define CONTROLLER_HPP

#include <cmath>
#include <utility>
#include <Eigen/Dense>
#include "Obstacle.hpp"
#include "CartesianFrenetConverter.hpp"

using namespace std;

class Controller {
    private:
        int N;
        double dt;
        double v;

        double w_d_mule;
        double w_d_hitch;
        double w_d_trailer;
        double w_d_dot_trailer;
        double w_phi;

        double w_omega_mule;
        double w_omega_rate_mule;

        double omega_min;
        double omega_max;
        double prev_omega;
        double max_delta_omega;

        bool use_exp_smooth;
        double smooth_alpha;

        bool is_reverse;

    public:
        Controller(int N = 12, double dt = 0.1, double v = 1.0,
            double w_d_mule = 1200.0, double w_d_hitch = 0.0, double w_d_trailer = 0.0,
            double w_d_dot_trailer = 200.0, double w_phi = 800.0,
            double w_omega_mule = 0.01, double w_omega_rate_mule = 1.0,
            pair<double, double> omega_bounds = {-2.0, 2.0}, double max_delta_omega = 0.5,
            bool use_exp_smooth = false, double smooth_alpha = 0.3, bool is_reverse = false);
        ~Controller();
        Eigen::Vector4d step_dynamics_forward(const Eigen::Vector4d& state, double omega, double L, double D) const;
        Eigen::Vector4d step_dynamics_reverse(const Eigen::Vector4d& state, double omega, double L, double D) const;
};

#endif
