#!/Users/subhodeep/venv/bin/python

import numpy as np
from Obstacle import Obstacle
from Controller import Controller
from ParametricFunctions import ParametricFunctions
from CartesianFrenetConverter import CartesianFrenetConverter
from PathPlanningFrenetFrame import PathPlanningFrenetFrame

if __name__ == "__main__":
    #========================================================================#
    #============= FORWARD FIGURE EIGHT WITH OBSTACLE AVOIDANCE =============#
    # direction toggle
    is_reverse = False

    path_func = lambda t: ParametricFunctions.figure_eight(t, a=10.0, v=1.0)
    controller = Controller(N=20, dt=0.1, v=1.0, is_reverse=is_reverse)
    sim = PathPlanningFrenetFrame(path_func, (0, 0, np.pi/4), controller, L=1.5, D=2.0, dt=0.1, T=80)

    # circular obstacle at (x,y) with radius r
    sim.add_obstacle(Obstacle(10.0, 0.0, 0.2),
                     clearance=1.0, influence_margin=1.0,
                     detour_span=5.0, pre_start=4.0, max_amplitude=0.5)
    sim.add_obstacle(Obstacle(-7.0, 5.0, 0.2),
                     clearance=1.0, influence_margin=1.0,
                     detour_span=5.0, pre_start=4.0, max_amplitude=0.5)

    sim.simulate()
    sim.animate()

    #========================================================================#
    #============ REVERSE FIGURE EIGHT WITHOUT OBSTACLE AVOIDANCE ===========#
    # direction toggle
    # is_reverse = True

    # path_func = lambda t: ParametricFunctions.figure_eight(t, a=10.0, v=1.0)
    # controller = Controller(N=20, dt=0.1, v=1.0, is_reverse=is_reverse)

    # x_start, y_start, theta_start, kappa_start = path_func(0.01)
    # start_pose = (x_start, y_start, CartesianFrenetConverter.normalize_angle(theta_start + np.pi))

    # sim = PathPlanningFrenetFrame(path_func, start_pose, controller, L=1.5, D=2.0, dt=0.1, T=80, is_reverse=is_reverse)

    # sim.simulate()
    # sim.animate()
