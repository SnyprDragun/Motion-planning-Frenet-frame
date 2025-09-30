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
    double dt = 0.1;                /*        sampling time [s]         */    
    int N = 20;                     /*        prediction horizon        */

    /* control limits */
    double v_min = -0.8;            /*    min linear velocity [m/s]     */
    double v_max = 1.0;             /*    max linear velocity [m/s]     */
    double w_min = -2.5;            /*   min angular velocity [rad/s]   */
    double w_max = 2.5;             /*   max angular velocity [rad/s]   */
    double dv_max = 0.5;            /*      max |Δv| per step [m/s]     */
    double dw_max = 1.0;            /*     max |Δω| per step [rad/s]    */

    /* costs */
    double q_xy = 2.0;              /*  weight on (x,y) position error  */
    double q_th = 0.3;              /*     weight on heading error      */
    double r_v = 0.02;              /*       weight on absolute v       */
    double r_w = 0.02;              /*       weight on absolute ω       */
    double s_v = 0.5;               /*          weight on |Δv|          */
    double s_w = 0.3;               /*          weight on |Δω|          */
    double qT_xy = 6.0;             /*          terminal (x,y)          */
    double qT_th = 1.0;             /*         terminal heading         */

    MPCParams(double dt_ = 0.1, int N_ = 15) : dt(dt_), N(N_) {}
};

class Controller {
    private:
        MatrixXd rollout(const VectorXd&, const VectorXd&, const VectorXd&);
        VectorXd warm_start(const VectorXd&, const VectorXd&);
        double normalize_angle(double);
        static double objective_wrapper(const vector<double>&, vector<double>&, void*);
        struct ObjData {
            Controller* mpc;
            VectorXd x_now;
            VectorXd x_target;
        };

    public:
        MPCParams p;
        VectorXd u_prev_seq;

        Controller(const MPCParams& params = MPCParams());
        vector<float> control(const PathPoint&, const PathPoint&);
};


#endif
