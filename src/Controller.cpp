/*
    Script for controller class
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#include "Controller.hpp"

Controller::Controller(const MPCParams& params) : p(params) {
    u_prev_seq = VectorXd::Zero(2 * p.N);
}

vector<float> Controller::control(const PathPoint& x_now, const PathPoint& x_target) {
    VectorXd x(3);
    x << x_now.x, x_now.y, x_now.theta;
    VectorXd x_ref(3);
    x_ref << x_target.x, x_target.y, x_target.theta;

    // Warm start
    VectorXd u0 = warm_start();

    // NLopt setup
    nlopt::opt opt(nlopt::LD_SLSQP, 2 * p.N);
    vector<double> lb(2*p.N), ub(2*p.N);
    for (int i=0;i<p.N;i++) { lb[i]=p.v_min; ub[i]=p.v_max; }
    for (int i=p.N;i<2*p.N;i++) { lb[i]=p.w_min; ub[i]=p.w_max; }
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);

    ObjData objdata{this, x, x_ref};
    opt.set_min_objective(objective_wrapper, &objdata);
    opt.set_maxeval(100);
    opt.set_ftol_rel(1e-4);

    vector<double> u_vec(u0.data(), u0.data()+u0.size());
    double minf;
    nlopt::result result = nlopt::FAILURE;
    try {
        result = opt.optimize(u_vec, minf);
    } catch(...) {
        u_vec.assign(u0.data(), u0.data()+u0.size()); // fallback
    }

    u_prev_seq = Map<VectorXd>(u_vec.data(), u_vec.size());

    auto clamp = [](float val, float min_val, float max_val) {
        return max(min_val, min(val, max_val));
    };
    float v_cmd = clamp(u_prev_seq(0), p.v_min, p.v_max);
    float w_cmd = clamp(u_prev_seq(p.N), p.w_min, p.w_max);
    cout << "Control commands: v = " << v_cmd << ", w = " << w_cmd << endl;
    return {v_cmd, w_cmd};
}

MatrixXd Controller::rollout(const VectorXd& x0, const VectorXd& v_seq, const VectorXd& w_seq) {
    MatrixXd traj(p.N+1, 3);
    VectorXd x = x0;
    traj.row(0) = x.transpose();
    for (int k=0;k<p.N;k++) {
        double v = v_seq(k);
        double w = w_seq(k);
        x(0) += p.dt * v * cos(x(2));
        x(1) += p.dt * v * sin(x(2));
        x(2) = normalize_angle(x(2) + p.dt*w);
        traj.row(k+1) = x.transpose();
    }
    return traj;
}

VectorXd Controller::warm_start() {
    VectorXd v_seq = u_prev_seq.head(p.N);
    VectorXd w_seq = u_prev_seq.tail(p.N);
    VectorXd v_ws(p.N), w_ws(p.N);
    v_ws.head(p.N-1) = v_seq.tail(p.N-1);
    v_ws(p.N-1) = v_seq(p.N-1);
    w_ws.head(p.N-1) = w_seq.tail(p.N-1);
    w_ws(p.N-1) = w_seq(p.N-1);
    VectorXd u_ws(2*p.N);
    u_ws << v_ws, w_ws;
    return u_ws;
}

double Controller::normalize_angle(double th) {
    while(th > M_PI) th -= 2*M_PI;
    while(th < -M_PI) th += 2*M_PI;
    return th;
}

double Controller::objective_wrapper(const vector<double>& u, vector<double>& grad, void* data) {
    ObjData* d = static_cast<ObjData*>(data);
    VectorXd u_vec = Map<const VectorXd>(u.data(), u.size());
    VectorXd v_seq = u_vec.head(d->mpc->p.N);
    VectorXd w_seq = u_vec.tail(d->mpc->p.N);

    MatrixXd traj = d->mpc->rollout(d->x_now, v_seq, w_seq);

    // Stage costs
    VectorXd pos_err2(traj.rows());
    for (int k=0;k<traj.rows();k++)
        pos_err2(k) = pow(traj(k,0)-d->x_target(0),2) + pow(traj(k,1)-d->x_target(1),2);

    VectorXd e_th2(traj.rows());
    for (int k=0;k<traj.rows();k++)
        e_th2(k) = pow(d->mpc->normalize_angle(traj(k,2)-d->x_target(2)),2);

    VectorXd dv = v_seq;
    dv.tail(dv.size()-1) -= v_seq.head(dv.size()-1);
    dv(0) = 0.0;
    VectorXd dw = w_seq;
    dw.tail(dw.size()-1) -= w_seq.head(dw.size()-1);
    dw(0) = 0.0;

    double J_stage = 0.0;
    for (int k=0;k<d->mpc->p.N;k++) {
        J_stage += d->mpc->p.q_xy * pos_err2(k) + d->mpc->p.q_th * e_th2(k)
                 + d->mpc->p.r_v * pow(v_seq(k),2) + d->mpc->p.r_w * pow(w_seq(k),2)
                 + d->mpc->p.s_v * abs(dv(k)) + d->mpc->p.s_w * abs(dw(k));
    }

    double J_term = d->mpc->p.qT_xy * pos_err2(d->mpc->p.N) + d->mpc->p.qT_th * e_th2(d->mpc->p.N);
    return J_stage + J_term;
}
