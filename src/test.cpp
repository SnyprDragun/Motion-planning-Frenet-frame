// #include <iostream>
// #include <fstream>
// #include <vector>
// #include <array>
// #include <cmath>
// #include <algorithm>

// using namespace std;
// using Point = array<double,2>;

// // -------- Utility Functions --------
// void write_csv(const string& filename, const vector<Point>& data) {
//     ofstream file(filename);
//     file << "x,y\n";
//     for (auto& p : data) file << p[0] << "," << p[1] << "\n";
//     file.close();
// }

// double norm(const Point& p) {
//     return sqrt(p[0]*p[0] + p[1]*p[1]);
// }

// Point normalize(const Point& p) {
//     double n = norm(p);
//     if (n < 1e-6) return {1.0, 0.0};
//     return {p[0]/n, p[1]/n};
// }

// Point operator+(const Point& a, const Point& b) { return {a[0]+b[0], a[1]+b[1]}; }
// Point operator-(const Point& a, const Point& b) { return {a[0]-b[0], a[1]-b[1]}; }
// Point operator*(double s, const Point& p) { return {s*p[0], s*p[1]}; }
// Point operator*(const Point& p, double s) { return {s*p[0], s*p[1]}; }

// // -------- Figure Eight Path --------
// Point figure_eight(double t) {
//     return {10.0 * sin(t), 10.0 * sin(t) * cos(t)};
// }

// vector<Point> generate_path(int N = 600, double tmax = 2*M_PI) {
//     vector<Point> path;
//     path.reserve(N);
//     for (int i = 0; i < N; i++) {
//         double t = tmax * i / (N-1);
//         path.push_back(figure_eight(t));
//     }
//     return path;
// }

// // -------- Simulation without Controller --------
// void simulate_three_point_body(const vector<Point>& path, vector<Point>& mule_traj,
//                                vector<Point>& joint_traj, vector<Point>& trailer_traj,
//                                double L=1.5, double joint_flex=0.3) {
//     Point trailer_dir = {-1.0, 0.0};
//     for (size_t i = 0; i < path.size(); i++) {
//         Point mule = path[i];
//         Point prev = (i > 0) ? path[i-1] : mule;
//         Point mule_dir = normalize(mule - prev);

//         // Flexible joint rotation
//         double angle_offset = joint_flex * sin(i * 0.05);
//         double c = cos(angle_offset), s = sin(angle_offset);
//         Point joint_dir = {-c*mule_dir[0] + s*mule_dir[1],
//                            -s*mule_dir[0] - c*mule_dir[1]};
//         Point joint = mule + L * joint_dir;

//         trailer_dir = normalize(0.9 * trailer_dir + 0.1 * joint_dir);
//         Point trailer = joint + L * trailer_dir;

//         mule_traj.push_back(mule);
//         joint_traj.push_back(joint);
//         trailer_traj.push_back(trailer);
//     }
// }

// // -------- Simple MPC-like Controller --------
// void mpc_controller(const vector<Point>& path, vector<Point>& mule_traj,
//                     vector<Point>& joint_traj, vector<Point>& trailer_traj,
//                     double L=1.5, double gain=0.25, double smooth=0.8) {
//     Point mule = path[0], joint = mule - Point{L, 0.0}, trailer = mule - Point{2*L, 0.0};
//     Point trailer_dir = {-1.0, 0.0}, mule_vel = {0.0, 0.0};

//     for (size_t i = 0; i < path.size(); i++) {
//         Point target = path[min(i+5, path.size()-1)];
//         Point desired_dir = normalize(target - mule);
//         mule_vel = smooth*mule_vel + (1.0 - smooth)*(gain * desired_dir);
//         mule = mule + mule_vel;

//         Point dir_joint = normalize(mule - joint);
//         joint = mule - L * dir_joint;

//         trailer_dir = normalize(0.85 * trailer_dir + 0.15 * dir_joint);
//         trailer = joint - L * trailer_dir;

//         mule_traj.push_back(mule);
//         joint_traj.push_back(joint);
//         trailer_traj.push_back(trailer);
//     }
// }

// // -------- Compute Error --------
// vector<double> compute_error(const vector<Point>& trailer_traj, const vector<Point>& path) {
//     vector<double> errors;
//     errors.reserve(trailer_traj.size());
//     for (auto& p : trailer_traj) {
//         double min_dist = 1e9;
//         for (auto& q : path) {
//             double dx = p[0] - q[0], dy = p[1] - q[1];
//             min_dist = min(min_dist, sqrt(dx*dx + dy*dy));
//         }
//         errors.push_back(min_dist);
//     }
//     return errors;
// }

// void write_error_csv(const string& filename, const vector<double>& err) {
//     ofstream file(filename);
//     file << "error\n";
//     for (auto e : err) file << e << "\n";
//     file.close();
// }

// // -------- Main --------
// int main() {
//     vector<Point> path = generate_path();

//     // No Controller
//     vector<Point> mule_nc, joint_nc, trailer_nc;
//     simulate_three_point_body(path, mule_nc, joint_nc, trailer_nc);
//     write_csv("mule_nc.csv", mule_nc);
//     write_csv("joint_nc.csv", joint_nc);
//     write_csv("trailer_nc.csv", trailer_nc);

//     // With Controller
//     vector<Point> mule_mpc, joint_mpc, trailer_mpc;
//     mpc_controller(path, mule_mpc, joint_mpc, trailer_mpc);
//     write_csv("mule_mpc.csv", mule_mpc);
//     write_csv("joint_mpc.csv", joint_mpc);
//     write_csv("trailer_mpc.csv", trailer_mpc);

//     // Error comparison
//     vector<double> err_nc = compute_error(trailer_nc, path);
//     vector<double> err_mpc = compute_error(trailer_mpc, path);
//     write_error_csv("err_nc.csv", err_nc);
//     write_error_csv("err_mpc.csv", err_mpc);

//     cout << "[INFO] Simulation complete. CSV files generated for visualization.\n";
//     return 0;
// }



// ==============================================================================
// ==============================================================================
// ==============================================================================
// ==============================================================================

// #include <iostream>
// #include <fstream>
// #include <vector>
// #include <array>
// #include <cmath>
// #include <algorithm>

// using namespace std;
// using Point = array<double,2>;

// // ---------- Utility ----------
// void write_csv(const string& filename, const vector<Point>& data) {
//     ofstream file(filename);
//     file << "x,y\n";
//     for (auto& p : data) file << p[0] << "," << p[1] << "\n";
//     file.close();
// }

// double norm(const Point& p) { return sqrt(p[0]*p[0] + p[1]*p[1]); }

// Point normalize(const Point& p) {
//     double n = norm(p);
//     if (n < 1e-6) return {1.0, 0.0};
//     return {p[0]/n, p[1]/n};
// }

// Point operator+(const Point& a, const Point& b) { return {a[0]+b[0], a[1]+b[1]}; }
// Point operator-(const Point& a, const Point& b) { return {a[0]-b[0], a[1]-b[1]}; }
// Point operator*(double s, const Point& p) { return {s*p[0], s*p[1]}; }
// Point operator*(const Point& p, double s) { return {s*p[0], s*p[1]}; }

// // ---------- Path Generation ----------
// Point figure_eight(double t) { return {10.0 * sin(t), 10.0 * sin(t) * cos(t)}; }

// vector<Point> generate_path(int N = 600, double tmax = 2*M_PI) {
//     vector<Point> path;
//     path.reserve(N);
//     for (int i = 0; i < N; i++) {
//         double t = tmax * i / (N-1);
//         path.push_back(figure_eight(t));
//     }
//     return path;
// }

// // ---------- Obstacle Struct ----------
// struct Obstacle { Point pos; double radius; };

// // ---------- Helper ----------
// inline double dist(const Point& a, const Point& b) {
//     return sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]));
// }

// Point low_pass(Point prev, Point new_dir, double smooth=0.8) {
//     return smooth*prev + (1.0 - smooth)*new_dir;
// }

// // ---------- Forward Motion ----------
// void simulate_three_point_body(const vector<Point>& path, vector<Point>& mule_traj,
//                                vector<Point>& joint_traj, vector<Point>& trailer_traj,
//                                double L=1.5, double joint_flex=0.3) {
//     Point trailer_dir = {-1.0, 0.0}, mule_vel={0,0};
//     for (size_t i = 0; i < path.size(); i++) {
//         Point mule = path[i];
//         if (i>0) mule_vel = low_pass(mule_vel, mule - path[i-1]);
//         mule = (i>0) ? path[i-1] + mule_vel : mule;

//         Point prev = (i > 0) ? path[i-1] : mule;
//         Point mule_dir = normalize(mule - prev);

//         double angle_offset = joint_flex * sin(i * 0.05);
//         double c = cos(angle_offset), s = sin(angle_offset);
//         Point joint_dir = {-c*mule_dir[0] + s*mule_dir[1], -s*mule_dir[0] - c*mule_dir[1]};
//         Point joint = mule + L * joint_dir;

//         trailer_dir = normalize(0.9 * trailer_dir + 0.1 * joint_dir);
//         Point trailer = joint + L * trailer_dir;

//         mule_traj.push_back(mule);
//         joint_traj.push_back(joint);
//         trailer_traj.push_back(trailer);
//     }
// }

// // ---------- MPC Controller ----------
// void mpc_controller(const vector<Point>& path, vector<Point>& mule_traj,
//                     vector<Point>& joint_traj, vector<Point>& trailer_traj,
//                     double L=1.5, double gain=0.25, double smooth=0.8) {
//     Point mule = path[0], joint = mule - Point{L, 0}, trailer = mule - Point{2*L,0};
//     Point trailer_dir = {-1.0, 0.0}, mule_vel={0,0};
//     for (size_t i=0; i<path.size(); i++) {
//         Point target = path[min(i+5, path.size()-1)];
//         Point dir = normalize(target - mule);
//         mule_vel = low_pass(mule_vel, gain * dir);
//         mule = mule + mule_vel;

//         Point dir_joint = normalize(mule - joint);
//         joint = mule - L * dir_joint;

//         trailer_dir = normalize(low_pass(trailer_dir, dir_joint));
//         trailer = joint - L * trailer_dir;

//         mule_traj.push_back(mule);
//         joint_traj.push_back(joint);
//         trailer_traj.push_back(trailer);
//     }
// }

// // ---------- Reverse Tracking ----------
// void simulate_reverse_tracking(const vector<Point>& path,
//                                vector<Point>& mule_traj,
//                                vector<Point>& joint_traj,
//                                vector<Point>& trailer_traj,
//                                double L = 1.5,
//                                double smooth = 0.85,
//                                double gain = 0.2) 
// {
//     vector<Point> rev_path(path.rbegin(), path.rend());

//     Point trailer = rev_path[0];
//     Point dir_init = normalize(rev_path[1] - rev_path[0]);
//     Point joint = trailer + L * dir_init;
//     Point mule  = joint + L * dir_init;

//     Point mule_vel = {0.0, 0.0};

//     for (size_t i = 0; i < rev_path.size(); i++) {
//         // 1. Trailer follows the reversed path
//         Point target_trailer = rev_path[min(i + 3, rev_path.size() - 1)];
//         Point dir_trailer = normalize(target_trailer - trailer);
//         trailer = trailer + 0.15 * dir_trailer;

//         // 2. Joint stays rigidly L ahead of trailer
//         Point push_dir = normalize(trailer - joint);
//         joint = trailer - L * push_dir;

//         // 3. Mule target position is L ahead of joint along push_dir
//         Point desired_mule = joint - L * push_dir;

//         // 4. Apply smoothed velocity towards desired_mule
//         mule_vel = low_pass(mule_vel, gain * normalize(desired_mule - mule), smooth);
//         mule = mule + mule_vel;

//         // 5. Enforce rigid distance constraint for mule (exactly L from joint)
//         mule = joint - L * push_dir;

//         // 6. Save trajectories
//         mule_traj.push_back(mule);
//         joint_traj.push_back(joint);
//         trailer_traj.push_back(trailer);
//     }
// }

// // ---------- Arc Generation & Smoothing ----------
// vector<Point> generate_arc_around_obstacle(const Point& entry, const Point& exit, const Obstacle& obs, double clearance);
// vector<Point> smooth_path(const vector<Point>& path, double max_angle_deg = 15.0);

// // (Use your previously defined arc and smooth functions unchanged here)

// // ---------- Generate Avoiding Path ----------
// vector<Point> generate_smooth_avoiding_path(const vector<Point>& original_path, const vector<Obstacle>& obstacles, double clearance = 1.0);

// // ---------- MPC with Avoidance ----------
// void mpc_with_avoidance(const vector<Point>& original_path, const vector<Obstacle>& obstacles,
//                         vector<Point>& mule_traj, vector<Point>& joint_traj, vector<Point>& trailer_traj,
//                         double L=1.5, double gain=0.25, double smooth=0.8) 
// {
//     vector<Point> safe_path = generate_smooth_avoiding_path(original_path, obstacles, 1.0);

//     // ✅ Save the generated safe path
//     write_csv("safe_path.csv", safe_path);

//     // ✅ Track this path using the existing MPC controller
//     mpc_controller(safe_path, mule_traj, joint_traj, trailer_traj, L, gain, smooth);
// }

// // ---------- Error ----------
// vector<double> compute_error(const vector<Point>& trailer_traj, const vector<Point>& path) {
//     vector<double> err(trailer_traj.size());
//     for (size_t i=0;i<trailer_traj.size();i++) {
//         double min_d=1e9;
//         for (auto& q:path) min_d=min(min_d,norm(trailer_traj[i]-q));
//         err[i]=min_d;
//     }
//     return err;
// }

// void write_error_csv(const string& fname,const vector<double>& err){
//     ofstream f(fname); f<<"error\n";
//     for(auto e:err) f<<e<<"\n";
//     f.close();
// }

// // ---- Arc Generation Around an Obstacle ----
// vector<Point> generate_arc_around_obstacle(const Point& entry, const Point& exit, const Obstacle& obs, double clearance) {
//     vector<Point> arc_points;
//     double R = obs.radius + clearance;

//     double a1 = atan2(entry[1] - obs.pos[1], entry[0] - obs.pos[0]);
//     double a2 = atan2(exit[1] - obs.pos[1], exit[0] - obs.pos[0]);

//     if (fabs(a2 - a1) > M_PI) {
//         if (a1 > a2) a2 += 2 * M_PI; else a1 += 2 * M_PI;
//     }

//     int N = 12; // more points → smoother arc
//     for (int i = 0; i <= N; i++) {
//         double a = a1 + (a2 - a1) * (i / (double)N);
//         arc_points.push_back({obs.pos[0] + R * cos(a), obs.pos[1] + R * sin(a)});
//     }
//     return arc_points;
// }

// // ---- Path Smoothing (Corner Rounding) ----
// vector<Point> smooth_path(const vector<Point>& path, double max_angle_deg) {
//     vector<Point> smoothed;
//     double max_angle = max_angle_deg * M_PI / 180.0;

//     smoothed.push_back(path.front());
//     for (size_t i = 1; i < path.size() - 1; i++) {
//         Point prev = path[i - 1], curr = path[i], next = path[i + 1];
//         Point v1 = normalize(curr - prev), v2 = normalize(next - curr);
//         double dotp = v1[0] * v2[0] + v1[1] * v2[1];
//         double angle = acos(clamp(dotp, -1.0, 1.0));

//         if (angle > max_angle) {
//             Point mid = {(curr[0] + next[0]) / 2.0, (curr[1] + next[1]) / 2.0};
//             smoothed.push_back(mid);
//         }
//         smoothed.push_back(curr);
//     }
//     smoothed.push_back(path.back());
//     return smoothed;
// }

// // ---- Final Smoothening: Moving Average + Chaikin-like Curve Smoothing ----
// #include <numeric> // for std::iota

// // ---- Step 2.5: Improved Final Smoothening using Cubic B-Spline ----
// // ---- Step 2.5: Localized Corner Rounding (Improved) ----
// vector<Point> final_path_smoothing(const vector<Point>& path, double max_angle_deg = 10.0, double corner_radius = 0.5) {
//     if (path.size() < 3) return path;

//     vector<Point> smooth_path;
//     double max_angle = max_angle_deg * M_PI / 180.0;

//     smooth_path.push_back(path.front());

//     for (size_t i = 1; i < path.size() - 1; i++) {
//         Point prev = path[i - 1], curr = path[i], next = path[i + 1];
//         Point v1 = normalize(curr - prev), v2 = normalize(next - curr);
//         double dotp = v1[0] * v2[0] + v1[1] * v2[1];
//         double angle = acos(clamp(dotp, -1.0, 1.0));

//         if (angle > max_angle) {
//             // ✅ Insert two intermediate points to round the corner
//             Point p1 = curr - v1 * corner_radius;
//             Point p2 = curr + v2 * corner_radius;
//             smooth_path.push_back(p1);
//             smooth_path.push_back(p2);
//         } else {
//             smooth_path.push_back(curr);
//         }
//     }

//     smooth_path.push_back(path.back());
//     return smooth_path;
// }


// // ---- Main Avoidance Path Generator ----
// vector<Point> generate_smooth_avoiding_path(const vector<Point>& original_path,
//                                             const vector<Obstacle>& obstacles,
//                                             double clearance) 
// {
//     vector<Point> new_path;
//     size_t i = 0;

//     while (i < original_path.size()) {
//         bool handled = false;

//         for (const auto& obs : obstacles) {
//             if (dist(original_path[i], obs.pos) < obs.radius + clearance) {
//                 size_t entry_idx = (i > 0) ? i - 1 : i;
//                 while (i < original_path.size() && dist(original_path[i], obs.pos) < obs.radius + clearance) i++;
//                 size_t exit_idx = (i < original_path.size()) ? i : original_path.size() - 1;

//                 auto arc_pts = generate_arc_around_obstacle(original_path[entry_idx], original_path[exit_idx], obs, clearance);
//                 new_path.insert(new_path.end(), arc_pts.begin(), arc_pts.end());
//                 handled = true;
//                 break;
//             }
//         }
//         if (!handled) new_path.push_back(original_path[i++]);
//     }

//     return final_path_smoothing(smooth_path(new_path, 5.0), 2);
// }



// // ---------- MAIN ----------
// int main(){
//     vector<Point> path=generate_path();

//     // Forward (NC)
//     vector<Point> m_nc,j_nc,t_nc;
//     simulate_three_point_body(path,m_nc,j_nc,t_nc);
//     write_csv("mule_nc.csv",m_nc); write_csv("joint_nc.csv",j_nc); write_csv("trailer_nc.csv",t_nc);

//     // Forward (MPC)
//     vector<Point> m_mpc,j_mpc,t_mpc;
//     mpc_controller(path,m_mpc,j_mpc,t_mpc);
//     write_csv("mule_mpc.csv",m_mpc); write_csv("joint_mpc.csv",j_mpc); write_csv("trailer_mpc.csv",t_mpc);

//     // Reverse Tracking
//     vector<Point> m_rev,j_rev,t_rev;
//     simulate_reverse_tracking(path,m_rev,j_rev,t_rev);
//     write_csv("mule_reverse.csv",m_rev); write_csv("joint_reverse.csv",j_rev); write_csv("trailer_reverse.csv",t_rev);

//     // Obstacles
//     vector<Obstacle> obstacles = {{{10.0, 0.0}, 0.5}, {{-10.0, 0.0}, 0.8}, {{5.0, 4.0}, 0.8}};

//     // Forward with Obstacle Avoidance
//     vector<Point> m_obs, j_obs, t_obs;
//     mpc_with_avoidance(path, obstacles, m_obs, j_obs, t_obs);

//     // ✅ Write both naming conventions for Python compatibility
//     write_csv("mule_avoid.csv", m_obs);
//     write_csv("joint_avoid.csv", j_obs);
//     write_csv("trailer_avoid.csv", t_obs);
//     write_csv("mule_obstacles.csv", m_obs);
//     write_csv("joint_obstacles.csv", j_obs);
//     write_csv("trailer_obstacles.csv", t_obs);

//     // Errors
//     write_error_csv("err_nc.csv", compute_error(t_nc, path));
//     write_error_csv("err_mpc.csv", compute_error(t_mpc, path));

//     cout<<"[INFO] ✅ Simulation complete with all CSVs generated.\n";
// }




// ==============================================================================
// ==============================================================================
// ==============================================================================
// ==============================================================================





// void mpc_with_obstacles(const vector<Point>& path,
//                                  vector<Point>& mule_traj,
//                                  vector<Point>& joint_traj,
//                                  vector<Point>& trailer_traj,
//                                  const vector<Obstacle>& obstacles,
//                                  double L = 1.5, double gain = 0.25, double smooth = 0.8) 
// {
//     Point mule = path[0], joint = mule - Point{L, 0}, trailer = mule - Point{2 * L, 0};
//     Point trailer_dir = {-1.0, 0.0}, mule_vel = {0, 0};

//     for (size_t i = 0; i < path.size(); i++) {
//         // --- 1. Normal path-following direction ---
//         Point target = path[min(i + 5, path.size() - 1)];
//         Point dir_path = normalize(target - mule);

//         // --- 2. Strong obstacle avoidance (only mule) ---
//         Point avoid_vec = {0, 0};
//         double avoid_w = 0.0;
//         for (const auto& o : obstacles) {
//             Point diff = mule - o.pos;          // vector away from obstacle
//             double d = norm(diff);
//             double influence = o.radius + 4.0;  // detection radius
//             if (d < influence) {
//                 avoid_vec = normalize(diff);
//                 // weight increases rapidly as mule gets close
//                 avoid_w = max(avoid_w, min(1.0, (influence - d) / influence * 0.6));
//                 break;
//             }
//         }

//         // --- 3. Combine path-following and avoidance ---
//         Point combined_dir;
//         if (avoid_w > 0.0) {
//             combined_dir = normalize((1.0 - avoid_w) * dir_path + avoid_w * avoid_vec);
//         } else {
//             combined_dir = dir_path;
//         }

//         // --- 4. Mule update with low-pass filtering ---
//         mule_vel = low_pass(mule_vel, gain * combined_dir, smooth);
//         mule = mule + mule_vel;

//         // --- 5. Joint simply follows mule ---
//         Point dir_joint = normalize(mule - joint);
//         joint = mule - L * dir_joint;

//         // --- 6. Trailer simply follows joint ---
//         trailer_dir = normalize(low_pass(trailer_dir, dir_joint));
//         trailer = joint - L * trailer_dir;

//         // --- 7. Store trajectories ---
//         mule_traj.push_back(mule);
//         joint_traj.push_back(joint);
//         trailer_traj.push_back(trailer);
//     }
// }




// ==============================================================================
// ==============================================================================
// ==============================================================================
// ==============================================================================




#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
using namespace std;

using Point = array<double,2>;

// ---------- Utility ----------
void write_csv(const string& filename, const vector<Point>& data) {
    ofstream file(filename);
    file << "x,y\n";
    for (auto& p : data) file << p[0] << "," << p[1] << "\n";
    file.close();
}

double norm(const Point& p) { return sqrt(p[0]*p[0] + p[1]*p[1]); }
Point normalize(const Point& p) {
    double n = norm(p);
    if (n < 1e-6) return {1.0, 0.0};
    return {p[0]/n, p[1]/n};
}
Point operator+(const Point& a, const Point& b) { return {a[0]+b[0], a[1]+b[1]}; }
Point operator-(const Point& a, const Point& b) { return {a[0]-b[0], a[1]-b[1]}; }
Point operator*(double s, const Point& p) { return {s*p[0], s*p[1]}; }
Point operator*(const Point& p, double s) { return {s*p[0], s*p[1]}; }
inline double dist(const Point& a, const Point& b) { return norm(a-b); }

// ---------- Path Generation ----------
Point figure_eight(double t) { return {10.0 * sin(t), 10.0 * sin(t) * cos(t)}; }

vector<Point> generate_path(int N = 600, double tmax = 2*M_PI) {
    vector<Point> path;
    path.reserve(N);
    for (int i = 0; i < N; i++) {
        double t = tmax * i / (N-1);
        path.push_back(figure_eight(t));
    }
    return path;
}

// ---------- Obstacle ----------
struct Obstacle { Point pos; double radius; };

// ---------- Kinematics ----------
struct State {
    double x;      // mule x
    double y;      // mule y
    double theta;  // mule heading
    double phi;    // trailer angle
};

State propagate(const State& s, double v, double omega, double L, double D, double dt) {
    State ns = s;
    ns.x     += v * cos(s.theta) * dt;
    ns.y     += v * sin(s.theta) * dt;
    ns.theta += omega * dt;
    double phi_dot = (v * sin(s.theta - s.phi) - L * omega * cos(s.theta - s.phi)) / D;
    ns.phi   += phi_dot * dt;
    return ns;
}

Point hitch_pos(const State& s, double L) {
    return { s.x - L * cos(s.theta), s.y - L * sin(s.theta) };
}
Point trailer_pos(const State& s, double L, double D) {
    Point H = hitch_pos(s, L);
    return { H[0] - D * cos(s.phi), H[1] - D * sin(s.phi) };
}

// ---------- Simulation with Kinematics ----------
void simulate_with_kinematics(const vector<Point>& path, vector<Point>& mule_traj,
                              vector<Point>& joint_traj, vector<Point>& trailer_traj,
                              double L=1.5, double D=2.0, double v=0.1, double gain=2.0, double dt=0.1) {
    State s{path[0][0], path[0][1], 0.0, 0.0};
    for (size_t i=1; i<path.size(); i++) {
        Point target = path[i];
        Point dir = normalize(target - Point{s.x, s.y});
        double desired_theta = atan2(dir[1], dir[0]);
        double theta_err = atan2(sin(desired_theta - s.theta), cos(desired_theta - s.theta));
        double omega = gain * theta_err;
        s = propagate(s, v, omega, L, D, dt);
        mule_traj.push_back({s.x, s.y});
        joint_traj.push_back(hitch_pos(s, L));
        trailer_traj.push_back(trailer_pos(s, L, D));
    }
}

// ---------- MPC (same dynamics) ----------
void mpc_controller_kinematics(const vector<Point>& path, vector<Point>& mule_traj,
                               vector<Point>& joint_traj, vector<Point>& trailer_traj,
                               double L=1.5, double D=2.0, double v=0.1, double gain=3.0, double dt=0.1) {
    State s{path[0][0], path[0][1], 0.0, 0.0};
    for (size_t i=0; i<path.size(); i++) {
        Point target = path[min(i+5, path.size()-1)];
        Point dir = normalize(target - Point{s.x, s.y});
        double desired_theta = atan2(dir[1], dir[0]);
        double theta_err = atan2(sin(desired_theta - s.theta), cos(desired_theta - s.theta));
        double omega = gain * theta_err;
        s = propagate(s, v, omega, L, D, dt);
        mule_traj.push_back({s.x, s.y});
        joint_traj.push_back(hitch_pos(s, L));
        trailer_traj.push_back(trailer_pos(s, L, D));
    }
}

// ---------- Reverse Tracking ----------
void simulate_reverse_kinematics(const vector<Point>& path, vector<Point>& mule_traj,
                                 vector<Point>& joint_traj, vector<Point>& trailer_traj,
                                 double L=1.5, double D=2.0, double v=-0.1, double gain=2.0, double dt=0.1) {
    vector<Point> rev_path(path.rbegin(), path.rend());
    State s{rev_path[0][0], rev_path[0][1], M_PI, M_PI}; // start facing backwards
    for (size_t i=1; i<rev_path.size(); i++) {
        Point target = rev_path[i];
        Point dir = normalize(target - Point{s.x, s.y});
        double desired_theta = atan2(dir[1], dir[0]);
        double theta_err = atan2(sin(desired_theta - s.theta), cos(desired_theta - s.theta));
        double omega = gain * theta_err;
        s = propagate(s, v, omega, L, D, dt);
        mule_traj.push_back({s.x, s.y});
        joint_traj.push_back(hitch_pos(s, L));
        trailer_traj.push_back(trailer_pos(s, L, D));
    }
}

// ---------- Obstacle Avoidance (reuse your old path generator) ----------
vector<Point> generate_arc_around_obstacle(const Point& entry, const Point& exit, const Obstacle& obs, double clearance) {
    vector<Point> arc_points;
    double R = obs.radius + clearance;
    double a1 = atan2(entry[1] - obs.pos[1], entry[0] - obs.pos[0]);
    double a2 = atan2(exit[1] - obs.pos[1], exit[0] - obs.pos[0]);
    if (fabs(a2 - a1) > M_PI) { if (a1 > a2) a2 += 2*M_PI; else a1 += 2*M_PI; }
    int N = 12;
    for (int i=0; i<=N; i++) {
        double a = a1 + (a2 - a1) * (i/(double)N);
        arc_points.push_back({obs.pos[0] + R*cos(a), obs.pos[1] + R*sin(a)});
    }
    return arc_points;
}
vector<Point> smooth_path(const vector<Point>& path, double max_angle_deg=15.0) {
    vector<Point> smoothed;
    double max_angle = max_angle_deg * M_PI / 180.0;
    smoothed.push_back(path.front());
    for (size_t i=1; i<path.size()-1; i++) {
        Point prev=path[i-1], curr=path[i], next=path[i+1];
        Point v1=normalize(curr-prev), v2=normalize(next-curr);
        double dotp=v1[0]*v2[0]+v1[1]*v2[1];
        double angle=acos(clamp(dotp,-1.0,1.0));
        if(angle>max_angle) smoothed.push_back((curr+next)*0.5);
        smoothed.push_back(curr);
    }
    smoothed.push_back(path.back());
    return smoothed;
}
vector<Point> final_path_smoothing(const vector<Point>& path, double max_angle_deg=10.0, double corner_radius=0.5) {
    if(path.size()<3) return path;
    vector<Point> out; out.push_back(path.front());
    double max_angle=max_angle_deg*M_PI/180.0;
    for(size_t i=1;i<path.size()-1;i++){
        Point prev=path[i-1],curr=path[i],next=path[i+1];
        Point v1=normalize(curr-prev),v2=normalize(next-curr);
        double dotp=v1[0]*v2[0]+v1[1]*v2[1];
        double angle=acos(clamp(dotp,-1.0,1.0));
        if(angle>max_angle){ out.push_back(curr-v1*corner_radius); out.push_back(curr+v2*corner_radius); }
        else out.push_back(curr);
    }
    out.push_back(path.back());
    return out;
}
vector<Point> generate_smooth_avoiding_path(const vector<Point>& original_path, const vector<Obstacle>& obstacles, double clearance=1.0) {
    vector<Point> new_path; size_t i=0;
    while(i<original_path.size()){
        bool handled=false;
        for(const auto& obs:obstacles){
            if(dist(original_path[i],obs.pos)<obs.radius+clearance){
                size_t entry=(i>0)?i-1:i;
                while(i<original_path.size()&&dist(original_path[i],obs.pos)<obs.radius+clearance) i++;
                size_t exit=(i<original_path.size())?i:original_path.size()-1;
                auto arc=generate_arc_around_obstacle(original_path[entry],original_path[exit],obs,clearance);
                new_path.insert(new_path.end(),arc.begin(),arc.end()); handled=true; break;
            }
        }
        if(!handled) new_path.push_back(original_path[i++]);
    }
    return final_path_smoothing(smooth_path(new_path,5.0),2);
}

// ---------- Error ----------
vector<double> compute_error(const vector<Point>& trailer_traj, const vector<Point>& path) {
    vector<double> err(trailer_traj.size());
    for (size_t i=0;i<trailer_traj.size();i++) {
        double min_d=1e9;
        for (auto& q:path) min_d=min(min_d,norm(trailer_traj[i]-q));
        err[i]=min_d;
    }
    return err;
}
void write_error_csv(const string& fname,const vector<double>& err){
    ofstream f(fname); f<<"error\n"; for(auto e:err) f<<e<<"\n"; f.close();
}

// ---------- MAIN ----------
int main(){
    vector<Point> path=generate_path();
    write_csv("reference_path.csv", path);

    // Forward (NC)
    vector<Point> m_nc,j_nc,t_nc;
    simulate_with_kinematics(path,m_nc,j_nc,t_nc);
    write_csv("mule_nc.csv",m_nc); write_csv("joint_nc.csv",j_nc); write_csv("trailer_nc.csv",t_nc);

    // Forward (MPC)
    vector<Point> m_mpc,j_mpc,t_mpc;
    mpc_controller_kinematics(path,m_mpc,j_mpc,t_mpc);
    write_csv("mule_mpc.csv",m_mpc); write_csv("joint_mpc.csv",j_mpc); write_csv("trailer_mpc.csv",t_mpc);

    // Reverse Tracking
    vector<Point> m_rev,j_rev,t_rev;
    simulate_reverse_kinematics(path,m_rev,j_rev,t_rev);
    write_csv("mule_reverse.csv",m_rev); write_csv("joint_reverse.csv",j_rev); write_csv("trailer_reverse.csv",t_rev);

    // Obstacles
    vector<Obstacle> obstacles = {{{10.0, 0.0}, 0.5}, {{-10.0, 0.0}, 0.8}, {{5.0, 4.0}, 0.8}};

    // Forward with Obstacle Avoidance
    vector<Point> m_obs, j_obs, t_obs;
    vector<Point> safe_path = generate_smooth_avoiding_path(path, obstacles, 1.0);
    mpc_controller_kinematics(safe_path,m_obs,j_obs,t_obs);
    write_csv("safe_path.csv", safe_path);
    write_csv("mule_avoid.csv", m_obs); write_csv("joint_avoid.csv", j_obs); write_csv("trailer_avoid.csv", t_obs);

    // Errors
    write_error_csv("err_nc.csv", compute_error(t_nc, path));
    write_error_csv("err_mpc.csv", compute_error(t_mpc, path));

    cout<<"[INFO] ✅ Simulation complete with all CSVs generated.\n";
}
