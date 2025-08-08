#!/Users/subhodeep/venv/bin/python

import numpy as np
from Controller import Controller
from ParametricFunctions import ParametricFunctions
from PathPlanningFrenetFrame import PathPlanningFrenetFrame

if __name__ == "__main__":
    path_func = lambda t: ParametricFunctions.figure_eight(t, a=10.0, v=1.0)
    controller = Controller(N=10, dt=0.1, v=1.0)
    sim = PathPlanningFrenetFrame(path_func, (0, 0, np.pi/4), controller)
    sim.simulate()
    sim.animate()
