// mule_frenet_sim.cpp
// Single-file port of your Python Controller + PathPlanningFrenetFrame
// Dependencies: Eigen (for vectors/matrices)
// Compile: g++ -O2 mule_frenet_sim.cpp -I /path/to/eigen -std=c++17 -o mule_sim
//
// Outputs: states.csv, hitches.csv, trailers.csv, frenet.csv
// (Plot these with Python/matplotlib if you want animation.)

// to run:
// rm *.csv
// g++ -std=c++17 test.cpp -I/opt/homebrew/include/eigen3 -o test
// ./test

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <functional>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

struct Obstacle {
    double x;
    double y;
    double radius;
    Obstacle(double x_=0, double y_=0, double r_=0): x(x_), y(y_), radius(r_) {}
};

static double wrap_angle(double a) {
    double res = std::fmod(a + M_PI, 2.0 * M_PI);
    if (res < 0) res += 2.0 * M_PI;
    return res - M_PI;
}

class Controller {
public:
    int N;
    double dt;
    double v;

    double w_d_mule;
    double w_d_hitch;
    double w_d_trailer;

    double w_omega;
    double w_omega_rate;

    double omega_min, omega_max;
    double prev_omega;
    double max_delta_omega;

    bool use_exp_smooth;
    double smooth_alpha;

    // RNG for optimizer
    std::mt19937_64 rng;

    Controller(int N_=12, double dt_=0.1, double v_=1.0,
               double w_d_mule_=1200.0, double w_d_hitch_=0.0, double w_d_trailer_=0.0,
               double w_omega_=0.01, double w_omega_rate_=1.0,
               std::pair<double,double> omega_bounds={-2.0,2.0}, double max_delta_omega_=0.5,
               bool use_exp_smooth_=false, double smooth_alpha_=0.3)
        : N(N_), dt(dt_), v(v_), w_d_mule(w_d_mule_), w_d_hitch(w_d_hitch_), w_d_trailer(w_d_trailer_),
          w_omega(w_omega_), w_omega_rate(w_omega_rate_), prev_omega(0.0),
          max_delta_omega(max_delta_omega_), use_exp_smooth(use_exp_smooth_), smooth_alpha(smooth_alpha_)
    {
        omega_min = omega_bounds.first;
        omega_max = omega_bounds.second;
        std::random_device rd;
        rng.seed(rd());
    }

    // step dynamics of mule/hitch/trailer kinematics (state = [x,y,theta,phi])
    Eigen::Vector4d step_dynamics(const Eigen::Vector4d& state, double omega, double L, double D) {
        double x = state[0], y = state[1], theta = state[2], phi = state[3];
        double dx = v * std::cos(theta);
        double dy = v * std::sin(theta);
        double dtheta = omega;
        double dphi = (v * std::sin(theta - phi) - L * omega * std::cos(theta - phi)) / D;

        double theta_new = wrap_angle(theta + dtheta * dt);
        double phi_new = wrap_angle(phi + dphi * dt);

        Eigen::Vector4d out;
        out << x + dx * dt, y + dy * dt, theta_new, phi_new;
        return out;
    }

    double lateral_deviation(double x, double y, double theta_ref, double rx, double ry) {
        return std::sin(theta_ref) * (x - rx) - std::cos(theta_ref) * (y - ry);
    }

    void compute_hitch_trailer_pose(const Eigen::Vector3d& mule_pose, double phi, double L, double D,
                                    Eigen::Vector3d& hitch_pose, Eigen::Vector3d& trailer_pose) {
        double x = mule_pose[0], y = mule_pose[1], theta = mule_pose[2];
        double x_hitch = x - L * std::cos(theta);
        double y_hitch = y - L * std::sin(theta);
        double theta_hitch = wrap_angle(theta + phi);
        double x_trailer = x_hitch - D * std::cos(theta_hitch);
        double y_trailer = y_hitch - D * std::sin(theta_hitch);
        double theta_trailer = theta_hitch;
        hitch_pose << x_hitch, y_hitch, theta_hitch;
        trailer_pose << x_trailer, y_trailer, theta_trailer;
    }

    // cost for a whole omega sequence
    double cost(const std::vector<double>& omega_seq, Eigen::Vector4d state0, double t0,
                std::function<Eigen::Vector3d(double)> path_func, double L, double D) {
        double J = 0.0;
        double prev_omega_local = prev_omega;

        for (int i = 0; i < N; ++i) {
            double omega = omega_seq[i];
            state0 = step_dynamics(state0, omega, L, D);
            double x = state0[0], y = state0[1], theta = state0[2], phi = state0[3];

            double t_pred = t0 + (i + 1) * dt;

            Eigen::Vector3d ref = path_func(t_pred); // rx, ry, theta_ref
            Eigen::Vector3d ref_hitch = path_func(std::max(0.0, t_pred - L / v));
            Eigen::Vector3d ref_trailer = path_func(std::max(0.0, t_pred - (L + D) / v));

            double d_mule = lateral_deviation(x, y, ref[2], ref[0], ref[1]);

            Eigen::Vector3d hitch_pose(0,0,0), trailer_pose(0,0,0);
            compute_hitch_trailer_pose(Eigen::Vector3d(x,y,theta), phi, L, D, hitch_pose, trailer_pose);

            double d_hitch = lateral_deviation(hitch_pose[0], hitch_pose[1], ref_hitch[2], ref_hitch[0], ref_hitch[1]);
            double d_trailer = lateral_deviation(trailer_pose[0], trailer_pose[1], ref_trailer[2], ref_trailer[0], ref_trailer[1]);

            J += w_d_mule * (d_mule * d_mule);
            J += w_d_hitch * (d_hitch * d_hitch);
            J += w_d_trailer * (d_trailer * d_trailer);

            J += w_omega * (omega * omega);
            J += w_omega_rate * ((omega - prev_omega_local) * (omega - prev_omega_local));

            double heading_err = std::atan2(std::sin(theta - ref[2]), std::cos(theta - ref[2]));
            J += 200.0 * (heading_err * heading_err);

            prev_omega_local = omega;
        }
        return J;
    }

    // simple shooting/hill-climb optimizer returning first command
    double mpc(double x, double y, double theta, double phi,
               std::function<Eigen::Vector3d(double)> path_func, double t0, double L, double D) {
        // initial guess: repeated prev_omega
        std::vector<double> best_seq(N, prev_omega);
        Eigen::Vector4d state0; state0 << x, y, theta, phi;
        double best_cost = cost(best_seq, state0, t0, path_func, L, D);

        std::normal_distribution<double> ndist(0.0, 0.5 * (omega_max - omega_min));
        const int iters = 120;  // tuning: more iter => better opt but slower

        for (int it = 0; it < iters; ++it) {
            // propose by perturbing some elements
            std::vector<double> cand = best_seq;
            int idxs = 1 + (rng() % std::max(1, N/3));
            for (int k=0; k<idxs; ++k) {
                int j = rng() % N;
                double perturb = ndist(rng);
                cand[j] = std::clamp(cand[j] + perturb, omega_min, omega_max);
            }

            double c = cost(cand, state0, t0, path_func, L, D);
            if (c < best_cost) {
                best_cost = c;
                best_seq = cand;
            }
        }

        double omega_raw = best_seq[0];
        double omega_out = omega_raw;
        if (use_exp_smooth) {
            omega_out = smooth_alpha * omega_raw + (1.0 - smooth_alpha) * prev_omega;
        }
        // clamp by max delta
        omega_out = std::clamp(omega_out, prev_omega - max_delta_omega, prev_omega + max_delta_omega);
        // final clip to global bounds
        omega_out = std::clamp(omega_out, omega_min, omega_max);
        prev_omega = omega_out;
        return omega_out;
    }
};

// small helper to map frenet detour d(s) into cartesian (rx,ry,rtheta) -> (x,y,theta_det).
// We'll map by shifting along path normal: n = (-sin(theta_ref), cos(theta_ref))
// theta_det is approximated as rtheta + atan(d_s) (small-angle approx)
Eigen::Vector3d frenet_to_cartesian_simple(double rx, double ry, double rtheta, double d, double d_s) {
    double nx = -std::sin(rtheta);
    double ny =  std::cos(rtheta);
    double x = rx + d * nx;
    double y = ry + d * ny;
    double theta_det = wrap_angle(rtheta + std::atan2(d_s, 1.0));
    return Eigen::Vector3d(x, y, theta_det);
}

class PathPlanningFrenetFrame {
public:
    std::function<Eigen::Vector3d(double)> path_func; // returns rx,ry,rtheta given time
    Controller& controller;
    double L, D, dt, T;
    int steps;

    double x, y, theta, phi; // mule state
    double v;

    // logging
    std::vector<Eigen::Vector2d> states;
    std::vector<Eigen::Vector2d> hitches;
    std::vector<Eigen::Vector2d> trailers;
    std::vector<Eigen::Vector2d> frenet; // s,d for mule

    std::vector<Obstacle> obstacles;
    double clearance;
    double influence_margin;
    double detour_span;
    double pre_start;
    double max_amplitude;

    // detour variables
    bool detour_active;
    double detour_s0, detour_s1, detour_s_center, detour_A;
    Obstacle detour_obs;

    double sim_time;
    bool debug;

    PathPlanningFrenetFrame(std::function<Eigen::Vector3d(double)> path_func_,
                            const Eigen::Vector3d& start,
                            Controller& controller_,
                            double L_=1.5, double D_=2.0, double dt_=0.1, double T_=100.0)
      : path_func(path_func_), controller(controller_), L(L_), D(D_), dt(dt_), T(T_), x(start[0]), y(start[1]),
        theta(start[2]), phi(start[2]), v(controller_.v),
        clearance(1.0), influence_margin(3.0), detour_span(6.0),
        pre_start(2.0), max_amplitude(2.5),
        detour_active(false), detour_s0(0), detour_s1(0), detour_s_center(0), detour_A(0),
        detour_obs(), sim_time(0.0), debug(false)
    {
        steps = static_cast<int>(std::round(T / dt));
    }

    void add_obstacle(const Obstacle& o, double clearance_=1.0, double influence_margin_=3.0,
                      double detour_span_=6.0, double pre_start_=2.0, double max_amplitude_=2.5) {
        obstacles.push_back(o);
        clearance = clearance_;
        influence_margin = influence_margin_;
        detour_span = detour_span_;
        pre_start = pre_start_;
        max_amplitude = max_amplitude_;
    }

    double _s_of_time(double t) const { return t * v; }
    Eigen::Vector3d _path_point(double t) const { return path_func(t); }
    double _distance_point_to_obs(double px, double py, const Obstacle& obs) const {
        return std::hypot(px - obs.x, py - obs.y);
    }

    // lookahead along path to see if any path point is within influence_range of obstacles
    std::pair<bool, std::pair<Obstacle,double>> _should_plan_detour(double current_t) {
        double lookahead = controller.N * controller.dt;
        int checks = std::max(1, static_cast<int>(std::floor(lookahead / dt)));
        for (int k=0; k<checks; ++k) {
            double t_check = current_t + k * dt;
            Eigen::Vector3d p = _path_point(t_check);
            for (const auto& obs: obstacles) {
                double inner = obs.radius + clearance;
                double outer = inner + influence_margin;
                double d = _distance_point_to_obs(p[0], p[1], obs);
                if (d <= outer) {
                    return {true, {obs, t_check}};
                }
            }
        }
        return {false, {Obstacle(), 0.0}};
    }

    // find best s_center (approx) nearest to obstacle
    std::pair<double,double> _find_s_center(double detect_t, const Obstacle& obs, double search_span=6.0) {
        double t0 = std::max(0.0, detect_t - search_span);
        double t1 = detect_t + search_span;
        int Nsearch = std::max(20, static_cast<int>((t1 - t0)/dt));
        double best_t = detect_t;
        double best_d = std::numeric_limits<double>::infinity();
        for (int i=0; i<=Nsearch; ++i) {
            double tt = t0 + (t1 - t0) * (static_cast<double>(i) / Nsearch);
            Eigen::Vector3d p = _path_point(tt);
            double d = _distance_point_to_obs(p[0], p[1], obs);
            if (d < best_d) {
                best_d = d;
                best_t = tt;
            }
        }
        double s_center = _s_of_time(best_t);
        if (debug) std::cout << "[center] best_t=" << best_t << " best_d=" << best_d << " s_center=" << s_center << "\n";
        return {s_center, best_d};
    }

    void _plan_quintic_detour(double detect_t, const Obstacle& obs) {
        auto [s_center, best_d] = _find_s_center(detect_t, obs);
        double half = 0.5 * detour_span;
        double s0 = std::max(0.0, s_center - half - pre_start);
        double s1 = s_center + half;

        double t_center = s_center / v;
        Eigen::Vector3d r = _path_point(t_center);
        double rx = r[0], ry = r[1], rtheta = r[2];

        double vx = obs.x - rx;
        double vy = obs.y - ry;
        double nx = -std::sin(rtheta);
        double ny = std::cos(rtheta);
        double dot = vx * nx + vy * ny;
        double sign = (dot > 0.0) ? -1.0 : 1.0;

        double trailer_margin = std::max(0.0, D * 0.1 + L * 0.05);
        double A_req = sign * (obs.radius + clearance + trailer_margin);
        double A = (std::abs(A_req) > max_amplitude) ? (std::copysign(max_amplitude, A_req)) : A_req;

        detour_active = true;
        detour_s0 = s0;
        detour_s1 = s1;
        detour_s_center = s_center;
        detour_A = A;
        detour_obs = obs;
        if (debug) std::cout << "[plan_quintic] s0=" << s0 << " s_center=" << s_center << " s1=" << s1 << " A=" << A << "\n";
    }

    void _maybe_deactivate_detour(double current_s) {
        if (!detour_active) return;
        if (current_s > detour_s1 + 1e-2) {
            if (debug) std::cout << "[done] current_s=" << current_s << " > s1=" << detour_s1 << " -> deactivate\n";
            detour_active = false;
            detour_s0 = detour_s1 = detour_s_center = detour_A = 0.0;
        }
    }

    // quintic phi and derivatives
    static double phi_s(double z) { return 10*z*z*z - 15*z*z*z*z + 6*z*z*z*z*z; }
    static double phi_d(double z) { return 30*z*z - 60*z*z*z + 30*z*z*z*z; }
    static double phi_dd(double z) { return 60*z - 180*z*z + 120*z*z*z; }

    // returns d, d_s, d_ss
    std::tuple<double,double,double> _quintic_bump_at_s(double s) const {
        if (!detour_active) return {0.0, 0.0, 0.0};
        if (s < detour_s0 - 1e-9 || s > detour_s1 + 1e-9) return {0.0, 0.0, 0.0};
        double Ls = std::max(1e-6, detour_s1 - detour_s0);
        double z = (s - detour_s0) / Ls;
        double A = detour_A;
        double p = phi_s(z);
        double pd = phi_d(z);
        double pdd = phi_dd(z);
        double d = A * p;
        double d_s = A * pd / Ls;
        double d_ss = A * pdd / (Ls * Ls);
        return {d, d_s, d_ss};
    }

    // map desired frenet bump into cartesian point
    Eigen::Vector3d avoidance_path_frenet(double t_pred) {
        Eigen::Vector3d base = _path_point(t_pred);
        double s = _s_of_time(t_pred);
        if (!detour_active) return base;
        if (s < detour_s0 - 1e-9 || s > detour_s1 + 1e-9) return base;
        auto [d, d_s, d_ss] = _quintic_bump_at_s(s);
        // map to cartesian by shifting along normal
        Eigen::Vector3d cart = frenet_to_cartesian_simple(base[0], base[1], base[2], d, d_s);
        return cart;
    }

    // Minimal cartesian_to_frenet used only for logging (approx): returns s_cond,d_cond
    std::pair<Eigen::Vector3d, Eigen::Vector3d> cartesian_to_frenet(double s_guess, double rx, double ry, double rtheta,
                                                                    double xq, double yq, double theta_q) {
        // approximate d as earlier lateral_deviation
        double d = std::sin(rtheta) * (xq - rx) - std::cos(rtheta) * (yq - ry);
        Eigen::Vector3d s_cond(s_guess, v, 0.0);
        Eigen::Vector3d d_cond(d, 0.0, 0.0);
        return {s_cond, d_cond};
    }

    void simulate() {
        for (int i = 0; i < steps; ++i) {
            double t = i * dt;
            double current_s = _s_of_time(t);

            if (!detour_active) {
                auto det = _should_plan_detour(t);
                if (det.first) {
                    const Obstacle& obs = det.second.first;
                    double detect_t = det.second.second;
                    _plan_quintic_detour(detect_t, obs);
                }
            }

            _maybe_deactivate_detour(current_s);

            auto path_func_local = [this](double t_pred) {
                return this->avoidance_path_frenet(t_pred);
            };

            double omega = controller.mpc(x, y, theta, phi, path_func_local, t, L, D);

            // integrate dynamics (explicit Euler)
            double dphi = (v * std::sin(theta - phi) - L * omega * std::cos(theta - phi)) / D;
            phi = wrap_angle(phi + dphi * dt);
            x += v * std::cos(theta) * dt;
            y += v * std::sin(theta) * dt;
            theta = wrap_angle(theta + omega * dt);

            sim_time += dt;

            Eigen::Vector3d hitch, trailer;
            Eigen::Vector3d mule_pose(x, y, theta);
            compute_hitch_trailer(mule_pose, phi, hitch, trailer);

            states.emplace_back(x, y);
            hitches.emplace_back(hitch[0], hitch[1]);
            trailers.emplace_back(trailer[0], trailer[1]);

            // log frenet using original parametric ref
            Eigen::Vector3d ref = path_func(t);
            double s_guess = current_s;
            auto [s_cond, d_cond] = cartesian_to_frenet(s_guess, ref[0], ref[1], ref[2], x, y, theta);
            frenet.emplace_back(s_cond[0], d_cond[0]);
        }
    }

    void compute_hitch_trailer(const Eigen::Vector3d& mule_pose, double phi_val,
                               Eigen::Vector3d& hitch_pose, Eigen::Vector3d& trailer_pose) const {
        double mx = mule_pose[0], my = mule_pose[1], mtheta = mule_pose[2];
        double xh = mx - L * std::cos(mtheta);
        double yh = my - L * std::sin(mtheta);
        double th = wrap_angle(mtheta + phi_val);
        double xt = xh - D * std::cos(th);
        double yt = yh - D * std::sin(th);
        hitch_pose << xh, yh, th;
        trailer_pose << xt, yt, th;
    }

    void save_csv() {
        {
            std::ofstream f("states.csv");
            f << "x,y\n";
            for (auto &p: states) f << p[0] << "," << p[1] << "\n";
        }
        {
            std::ofstream f("hitches.csv");
            f << "x,y\n";
            for (auto &p: hitches) f << p[0] << "," << p[1] << "\n";
        }
        {
            std::ofstream f("trailers.csv");
            f << "x,y\n";
            for (auto &p: trailers) f << p[0] << "," << p[1] << "\n";
        }
        {
            std::ofstream f("frenet.csv");
            f << "s,d\n";
            for (auto &p: frenet) f << p[0] << "," << p[1] << "\n";
        }
        std::cout << "Saved CSVs: states.csv, hitches.csv, trailers.csv, frenet.csv\n";
    }
};

///// Example parametric figure-eight function similar to your Python one
Eigen::Vector3d figure_eight(double t, double a=10.0, double v=1.0) {
    // A simple parametric figure-eight (lemniscate-like)
    // We'll parametrize with angle phi = (v / a) * t
    double phi = (v / a) * t;
    double x = a * std::sin(phi);
    double y = a * std::sin(phi) * std::cos(phi); // somewhat figure-eight-ish
    // approximate heading as derivative direction
    double dx = a * std::cos(phi) * (v / a);
    double dy = (a * (std::cos(2*phi) * (v / a))); // derivative approx
    double theta = std::atan2(dy, dx);
    return Eigen::Vector3d(x, y, theta);
}

int main() {
    // start pose and controller
    Eigen::Vector3d start(0.0, 0.0, M_PI/4.0);
    Controller controller(10, 0.1, 1.0);
    controller.w_d_mule = 1200.0;
    controller.w_omega = 0.01;
    controller.w_omega_rate = 1.0;

    // create sim
    auto path_fn = [](double t) { return figure_eight(t, 10.0, 1.0); };
    PathPlanningFrenetFrame sim(path_fn, start, controller, 1.5, 2.0, 0.1, 60.0);

    // add obstacles (two as in python)
    sim.add_obstacle(Obstacle(10.0, 0.0, 0.2), 1.0, 1.0, 5.0, 4.0, 0.5);
    sim.add_obstacle(Obstacle(-7.0, 5.0, 0.2), 1.0, 1.0, 5.0, 4.0, 0.5);

    // enable debug if desired
    // sim.debug = true;

    std::cout << "Simulating...\n";
    sim.simulate();
    sim.save_csv();
    std::cout << "Done\n";
    return 0;
}
