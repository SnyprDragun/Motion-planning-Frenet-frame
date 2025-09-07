#!/Users/subhodeep/venv/bin/python
'''class for control loop'''
import numpy as np
from tqdm import tqdm
from tabulate import tabulate
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class PathPlanningFrenetFrame:
    def __init__(self, robot, target_path, controller, T, dt):
        self.robot = robot
        self.target_path = target_path
        self.controller = controller
        self.T = T
        self.dt = dt

        self.current_states = []
        self.target_states = []
        self.control_actions = []

        self.hitches = []
        self.trailers = []

        self.current_states_frenet = []
        self.target_states_frenet = []
        self.control_actions_frenet = []

    def store_hitch_trailer(self):
        hitches = []
        trailers = []

        for i in range(self.robot.trailer_count):
            hitch, trailer = self.robot.calculate_kth_hitch_trailer_pose(i + 1)
            hitches.append(hitch)
            trailers.append(trailer)

        self.hitches.append(hitches)
        self.trailers.append(trailers)

    def control_loop(self):
        current_state = np.append(self.robot.mule_position, self.robot.mule_orientation)

        for i in tqdm(range(int(self.T / self.dt)), desc="Simulating", unit="iterations"):
            t = i * self.dt

            target_state = self.target_path.equation(t)
            control_action = self.controller.control(current_state, target_state)
            updated_state = self.robot.update_state(control_action, self.dt)
            current_state = updated_state

            self.store_hitch_trailer()
            self.current_states.append(current_state)
            self.target_states.append(target_state)
            self.control_actions.append(control_action)

    def diagnostics(self):
        r'''
        Displays control step in tabular format
        +---------------------+-----------------------+-------------------+---....
        |     Current Pose    |      Target Pose      |   Control Action  |
        +=====================+=======================+===================+===....
        | [$x_{mule}$, $y_{mule}$, $\theta_{mule}$] | [$x_{target}$, $y_{target}$, $\theta_{target}$] |       [$v$, $\omega$]      |
        +---------------------+-----------------------+-------------------+---....
        '''
        headers = ["Sl. no"] + ["Current Pose"] + ["Target Pose"] + ["Control Action"]
        rows = []
        for i in range(int(self.T / self.dt)):
            row = [
                i + 1,
                [f"{x:.4f}" for x in np.array(self.current_states[i]).flatten()],
                [f"{x:.4f}" for x in np.array(self.target_states[i]).flatten()],
                [f"{x:.4f}" for x in np.array(self.control_actions[i]).flatten()]
            ]
            rows.append(row)
        print(tabulate(rows, headers, tablefmt='grid'))

    def plot(self, interval=100):
        self.target_states = np.array(self.target_states)
        self.current_states = np.array(self.current_states)
        self.control_actions = np.array(self.control_actions)
        self.hitches = np.array(self.hitches)
        self.trailers = np.array(self.trailers)

        n_frames = int(self.T / self.dt)
        t = np.arange(n_frames)

        fig = plt.figure(figsize=(12, 6))
        gs = fig.add_gridspec(2, 2, width_ratios=[2, 1])
        ax_traj = fig.add_subplot(gs[:, 0])
        ax_v = fig.add_subplot(gs[0, 1])
        ax_w = fig.add_subplot(gs[1, 1])

        # ---- Trajectory plot ---- #
        ax_traj.set_title("Trajectory")
        ax_traj.set_xlabel("X")
        ax_traj.set_ylabel("Y")
        ax_traj.plot(self.target_states[:, 0], self.target_states[:, 1], "k--", label="Ideal")
        current_line, = ax_traj.plot([], [], "r-", label="Simulated")
        ax_traj.legend()
        ax_traj.set_aspect("equal", adjustable="datalim")

        # ---- Control plots ---- #
        ax_v.set_title("Linear Velocity v")
        ax_v.set_xlabel("Time")
        ax_v.set_ylabel("v")
        line_v, = ax_v.plot([], [], "b-")

        ax_w.set_title("Angular Velocity ω")
        ax_w.set_xlabel("Time")
        ax_w.set_ylabel("ω")
        line_w, = ax_w.plot([], [], "g-")

        # ---- Hitch & Trailer plots ---- #
        n_pairs = self.hitches.shape[1]
        mule_scats = [ax_traj.plot([], [], "ro", label=f"Mule")[0]]
        hitch_scats = [ax_traj.plot([], [], "bo", label=f"Hitch {i+1}" if i == 0 else "")[0] for i in range(n_pairs)]
        trailer_scats = [ax_traj.plot([], [], "go", label=f"Trailer {i+1}" if i == 0 else "")[0] for i in range(n_pairs)]

        # ---- Links ---- #
        n_links = 2 * n_pairs
        link_lines = [ax_traj.plot([], [], "k-")[0] for _ in range(n_links)]

        def update(frame):
            current_line.set_data(self.current_states[:frame, 0], self.current_states[:frame, 1])
            line_v.set_data(t[:frame], self.control_actions[:frame, 0])
            line_w.set_data(t[:frame], self.control_actions[:frame, 1])

            mx, my = self.current_states[frame-1, 0], self.current_states[frame-1, 1]
            mule_scats[0].set_data([mx], [my])

            for i in range(n_pairs):
                hx, hy = self.hitches[frame, i]
                tx, ty = self.trailers[frame, i]
                hitch_scats[i].set_data([hx], [hy])
                trailer_scats[i].set_data([tx], [ty])

            chain_x = [mx]
            chain_y = [my]
            for i in range(n_pairs):
                hx, hy = self.hitches[frame, i]
                tx, ty = self.trailers[frame, i]
                chain_x += [hx, tx]
                chain_y += [hy, ty]

            for j in range(len(chain_x) - 1):
                link_lines[j].set_data([chain_x[j], chain_x[j+1]], [chain_y[j], chain_y[j+1]])

            ax_traj.relim(); ax_traj.autoscale_view()
            ax_v.relim(); ax_v.autoscale_view()
            ax_w.relim(); ax_w.autoscale_view()

            return [current_line, line_v, line_w] + mule_scats + hitch_scats + trailer_scats + link_lines

        ani = animation.FuncAnimation(fig, update, frames=n_frames, interval=interval, blit=False, repeat=False)
        plt.show()

    def display_time(self, start, end):
        k = int(end - start)
        mins = (k // 60)
        if end - start < 1:
            secs = (((end - start) * 10000) // 100) / 100
        else:
            secs = k - (mins * 60)
        print(f"Time taken: {mins} minutes {secs} seconds")
