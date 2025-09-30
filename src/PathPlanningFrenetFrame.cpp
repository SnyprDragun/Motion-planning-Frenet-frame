/*
    Script for path planning and task solver class
    Author: Subhodeep Choudhury
    Source: https://github.com/SnyprDragun/Motion-planning-Frenet-frame
*/

#include "PathPlanningFrenetFrame.hpp"

ostream& operator<<(ostream& os, const PathPoint& pt) {
    os << pt.x << ", " << pt.y << ", " << pt.theta;
    return os;
}

PathPlanningFrenetFrame::PathPlanningFrenetFrame(RobotDynamics& robot, Path& target_path, Controller& controller, float T, float dt
) : robot(robot), target_path(target_path), controller(controller), T(T), dt(dt), current_states(), target_states(), 
control_actions(), hitches(), trailers(), current_states_frenet(), target_states_frenet(), control_actions_frenet() {}

PathPlanningFrenetFrame::~PathPlanningFrenetFrame(){}

void PathPlanningFrenetFrame::store_hitch_trailer() {
    vector<vector<float>> hitches_step;
    vector<vector<float>> trailers_step;

    for (int i = 0; i < this->robot.trailer_count; ++i) {
        auto [hitch, trailer] = this->robot.calculate_kth_hitch_trailer_pose(i + 1);
        hitches_step.push_back(hitch);
        trailers_step.push_back(trailer);
    }

    this->hitches.push_back(hitches_step);
    this->trailers.push_back(trailers_step);
}

void PathPlanningFrenetFrame::control_loop() {
    PathPoint current_state;
    current_state.x = this->robot.mule_position[0];
    current_state.y = this->robot.mule_position[1];
    current_state.theta = this->robot.mule_orientation;

    int n_iterations = static_cast<int>(T / dt);

    for (int i = 0; i < n_iterations; ++i) {
        float t = i * dt;

        PathPoint target_state = this->target_path.equation(t);
        vector<float> control_action = this->controller.control(current_state, target_state);
        PathPoint updated_state = this->robot.update_state(control_action, dt);
        current_state = updated_state;

        store_hitch_trailer();
        this->current_states.push_back(current_state);
        this->target_states.push_back(target_state);
        this->control_actions.push_back(control_action);
    }
}

void PathPlanningFrenetFrame::diagnostics() {
    int n_steps = static_cast<int>(T / dt);

    cout << setw(8) << "Sl.no"
              << setw(25) << "Current Pose"
              << setw(25) << "Target Pose"
              << setw(20) << "Control Action" << endl;

    cout << string(78, '-') << endl;

    for (int i = 0; i < n_steps; ++i) {
        cout << setw(8) << i + 1;

        cout << setw(25) << "[";
        cout << fixed << setprecision(4) << current_states[i];
        cout << "]";
        cout << setw(25) << "[";
        cout << fixed << setprecision(4) << target_states[i];
        cout << "]";
        cout << setw(20) << "[";

        for (size_t j = 0; j < control_actions[i].size(); ++j) {
            cout << fixed << setprecision(4) << control_actions[i][j];
            if (j != control_actions[i].size() - 1) cout << ", ";
        }
        cout << "]";

        cout << endl;
    }
}

void PathPlanningFrenetFrame::plot(int interval_ms) {
    int n_frames = static_cast<int>(T / dt);

    vector<double> target_x, target_y;
    for (auto& state : target_states) {
        target_x.push_back(state.x);
        target_y.push_back(state.y);
    }

    plt::figure_size(1200, 600);

    for (int frame = 0; frame < n_frames; ++frame) {
        vector<double> current_x, current_y;
        for (int i = 0; i <= frame; ++i) {
            current_x.push_back(current_states[i].x);
            current_y.push_back(current_states[i].y);
        }

        plt::clf();
        plt::plot(target_x, target_y, "k--"); 
        plt::plot(current_x, current_y, "r-"); 
        plt::plot({current_states[frame].x}, {current_states[frame].y}, "ro");

        int n_pairs = hitches[0].size();
        vector<double> link_x = {current_states[frame].x};
        vector<double> link_y = {current_states[frame].y};

        for (int i = 0; i < n_pairs; ++i) {
            double hx = hitches[frame][i][0];
            double hy = hitches[frame][i][1];
            double tx = trailers[frame][i][0];
            double ty = trailers[frame][i][1];

            plt::plot({hx}, {hy}, "bo");
            plt::plot({tx}, {ty}, "go");

            link_x.push_back(hx);
            link_y.push_back(hy);
            link_x.push_back(tx);
            link_y.push_back(ty);
        }

        plt::plot(link_x, link_y, "k-");
        plt::xlim(-12, 12);
        plt::ylim(-12, 12);
        plt::pause(interval_ms / 1000.0);
    }

    plt::show();
}

void PathPlanningFrenetFrame::display_time(float start, float end) {
    float elapsed = end - start;
    int k = static_cast<int>(elapsed);
    int mins = k / 60;
    float secs;

    if (elapsed < 1.0) {
        secs = floor(elapsed * 10000.0) / 100.0;
    } else {
        secs = k - (mins * 60);
    }

    cout << "Time taken: " << mins << " minutes " << secs << " seconds" << endl;
}
