#!/Users/subhodeep/venv/bin/python

import numpy as np
from Obstacle import Obstacle
from Controller import Controller
from ParametricFunctions import ParametricFunctions
from PathPlanningFrenetFrame import PathPlanningFrenetFrame

if __name__ == "__main__":
    path_func = lambda t: ParametricFunctions.figure_eight(t, a=10.0, v=1.0)
    controller = Controller(N=10, dt=0.1, v=1.0)
    sim = PathPlanningFrenetFrame(path_func, (0, 0, np.pi/4), controller, L=1.5, D=2.0, dt=0.1, T=70)

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
