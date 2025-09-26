/*
    Header for controller class
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#ifndef CONTROLLER_HPP
#define CONTROLLER_HPP

#pragma once
#include <Eigen/Dense>
#include <vector>
#include <nlopt.hpp>
#include <cmath>
#include <iostream>
#include <algorithm>

#include "PathPoint.hpp"

using namespace std;
using namespace Eigen;

struct MPCParams {
    double dt = 0.1;
    int N = 15;

    double v_min = -0.8;
    double v_max = 1.0;
    double w_min = -2.5;
    double w_max = 2.5;

    double dv_max = 0.5;
    double dw_max = 1.0;

    double q_xy = 2.0;
    double q_th = 0.3;
    double r_v = 0.02;
    double r_w = 0.02;
    double s_v = 0.5;
    double s_w = 0.3;
    double qT_xy = 6.0;
    double qT_th = 1.0;

    MPCParams(double dt_ = 0.1, int N_ = 15) : dt(dt_), N(N_) {}
};

class Controller {
public:
    MPCParams p;
    VectorXd u_prev_seq; // [v0..vN-1, w0..wN-1]

    Controller(const MPCParams& params = MPCParams());

    // Main control function: returns (v, w)
    vector<float> control(const PathPoint& x_now, const PathPoint& x_target);

private:
    MatrixXd rollout(const VectorXd& x0, const VectorXd& v_seq, const VectorXd& w_seq);
    VectorXd warm_start();
    double normalize_angle(double th);

    // Static wrapper for NLopt
    static double objective_wrapper(const vector<double>& u, vector<double>& grad, void* data);
    struct ObjData {
        Controller* mpc;
        VectorXd x_now;
        VectorXd x_target;
    };
};


#endif
