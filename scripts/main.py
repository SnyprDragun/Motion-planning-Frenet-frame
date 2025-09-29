#!/Users/subhodeep/venv/bin/python3
'''main script'''
import time
import numpy as np
from Path import Path
from Obstacle import Obstacle
from RobotDynamics import RobotDynamics
from Controller import MPCParams, UnicycleMPC
from PathPlanningFrenetFrame import PathPlanningFrenetFrame

if __name__ == "__main__":
    r'''
    < ================= Main ================= >

    Choose:
        > robot dynamics - no. of trailers, direction of motion
        > controller - mpc or other
        > path - straight line, circle, ellipse, figure-eight

    Specify:
        > initial position of last trailer (offset)
        > initial values of all angles ($\theta$, $\phi$)
        > lengths of joining rods ($L$, $D$)
    '''
    start = time.time()

    OFFSET = [10,-7]
    THETA = [np.pi/2, np.pi/2]
    PHI = [np.pi/2, np.pi/2]
    L = [1.5, 1.5]
    D = [2.0, 2.0]

    robot = RobotDynamics(OFFSET, L, D, THETA, PHI, trailer_count=2, direction=True)
    robot.diagnostics()
    mpc = UnicycleMPC(MPCParams(dt=0.1, N=20))
    circle = Path("circle")
    trajectory = PathPlanningFrenetFrame(robot=robot, target_path=circle, controller=mpc, T=65, dt=0.1)
    trajectory.control_loop()
    # trajectory.diagnostics()
    robot.diagnostics()
    trajectory.display_time(start, time.time())
    trajectory.plot()
