/*
    Demo: Using CartesianFrenetConverter to convert a parametric curve 
    into Frenet coordinates and generate data for plotting.

    Author: Subhodeep Choudhury
*/

// #include "CartesianFrenetConverter.hpp"
// #include "ParametricFunctions.hpp"
// #include <iostream>
// #include <fstream>

// // -------- Parametric Curve (Figure Eight) --------
// array<double,2> figure_eight(double t) {
//     double x = 10.0 * sin(t);
//     double y = 10.0 * sin(t) * cos(t);
//     return {x, y};
// }

// // -------- Compute heading angle --------
// double compute_heading(const array<double,2>& p1, const array<double,2>& p2) {
//     return atan2(p2[1] - p1[1], p2[0] - p1[0]);
// }

// // -------- Main Program --------
// int main() {
//     const int N = 400;
//     const double T_MAX = 2 * M_PI;

//     vector<array<double,2>> points;
//     points.reserve(N);

//     ParametricFunctions param = ParametricFunctions(10);

//     // Sample the parametric curve
//     for (int i = 0; i < N; i++) {
//         double t = T_MAX * i / (N - 1);
//         points.push_back(param.figure_eight(t));
//     }

//     // Compute arc length s
//     vector<double> s_vals(N, 0.0);
//     for (int i = 1; i < N; i++) {
//         double dx = points[i][0] - points[i-1][0];
//         double dy = points[i][1] - points[i-1][1];
//         s_vals[i] = s_vals[i-1] + sqrt(dx*dx + dy*dy);
//     }

//     // Store Frenet coordinates
//     vector<double> frenet_s, frenet_d;
//     frenet_s.reserve(N);
//     frenet_d.reserve(N);

//     // Convert to Frenet frame using your class
//     for (int i = 0; i < N; i++) {
//         double rs = s_vals[i];
//         double rtheta = (i < N-1) ? compute_heading(points[i], points[i+1])
//                                   : compute_heading(points[i-1], points[i]);

//         double rkappa = 0.0;   // assume path curvature is negligible
//         double rdkappa = 0.0;  // assume curvature rate is negligible

//         auto s_d_pair = CartesianFrenetConverter::cartesian_to_frenet(
//             rs, points[i][0], points[i][1], rtheta, rkappa, rdkappa,
//             points[i][0], points[i][1], 1.0, 0.0, rtheta, 0.0);

//         frenet_s.push_back(s_d_pair.first[0]);   // s(t)
//         frenet_d.push_back(s_d_pair.second[0]);  // d(s)
//     }

//     // -------- Write Data to CSV --------
//     ofstream file("curve_frenet_output.csv");
//     file << "x,y,s,d\n";
//     for (int i = 0; i < N; i++) {
//         file << points[i][0] << "," << points[i][1] << ","
//              << frenet_s[i] << "," << frenet_d[i] << "\n";
//     }
//     file.close();

//     cout << "Data saved to curve_frenet_output.csv. Use Python or gnuplot to plot.\n";
//     return 0;
// }


#include "PathPlanningFrenetFrame.hpp"

int main() {
    PathPlanningFrenetFrame planner;

    planner.convert_cartesian_path_to_frenet();
    planner.track_frenet_path_no_obstacles();
    planner.track_frenet_path_with_obstacles();
    planner.track_frenet_path_with_trailer_no_obstacles();
    planner.track_frenet_path_with_trailer_with_obstacles_reverse();

    return 0;
}
