/*
    Header for conversions between Cartesian and Frenet-Serret frame
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#ifndef PATH_PLANNING_FRENET_FRAME_HPP
#define PATH_PLANNING_FRENET_FRAME_HPP

#include <vector>
#include <cmath>
#include <array>

using namespace std;

class PathPlanningFrenetFrame {
    private:
        

    public:
        PathPlanningFrenetFrame();
        ~PathPlanningFrenetFrame();
        void convert_cartesian_path_to_frenet();
        void track_frenet_path_no_obstacles();
        void track_frenet_path_with_obstacles();
        void track_frenet_path_with_trailer_no_obstacles();
        void track_frenet_path_with_trailer_with_obstacles_reverse();
};

#endif
