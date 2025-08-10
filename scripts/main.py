#!/Users/subhodeep/venv/bin/python

import numpy as np
from Obstacle import Obstacle
from Controller import Controller
from ParametricFunctions import ParametricFunctions
from CartesianFrenetConverter import CartesianFrenetConverter
from PathPlanningFrenetFrame import PathPlanningFrenetFrame

if __name__ == "__main__":
    # direction toggle
    is_reverse = False

    path_func = lambda t: ParametricFunctions.figure_eight(t, a=10.0, v=1.0)
    x_start, y_start, theta_start, kappa_start = path_func(0.0)
    # For reverse tracking, the trailer must face the opposite direction of the path's tangent
    start_pose = (x_start, y_start, CartesianFrenetConverter.normalize_angle(theta_start + np.pi))

    controller = Controller(N=20, dt=0.1, v=1.0, is_reverse=is_reverse)
    sim = PathPlanningFrenetFrame(path_func, (0, 0, np.pi/4), controller, L=1.5, D=2.0, dt=0.1, T=80)

    # Add a circular obstacle at (x,y) radius r
    sim.add_obstacle(Obstacle(10.0, 0.0, 0.2),
                     clearance=1.0, influence_margin=1.0,
                     detour_span=5.0, pre_start=4.0, max_amplitude=0.5)
    sim.add_obstacle(Obstacle(-7.0, 5.0, 0.2),
                     clearance=1.0, influence_margin=1.0,
                     detour_span=5.0, pre_start=4.0, max_amplitude=0.5)

    # Enable debug prints if you need internal info
    # sim.debug = True

    sim.simulate()
    sim.animate()
