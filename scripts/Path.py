#!/Users/subhodeep/venv/bin/python
'''class to generate path points for motion'''
import numpy as np

class Path():
    '''
    Parametric functions to generate path points. Also handles obstcale avoidance. 
    Returns. safe and optimal path point for traversal
    '''
    def __init__(self, shape, R=10.0, a=12.0, b=8.0, slope=0.0, intercept=0.0, v=1.0):
        self.shape = shape
        self.R = R
        self.a = a
        self.b = b
        self.slope = slope
        self.intercept = intercept
        self.v = v
        self.obstacles = 0

    @staticmethod
    def circle(t, R=10.0, v=1.0):
        r'''
        Generates a circle path and its kinematic properties.
        Returns: $x$, $y$, $\theta$ (heading), $\kappa$ (curvature)
        • reaches back at starting point at $t = R * 2\pi$ (62.8)
        '''
        omega = v / R
        x = R * np.cos(omega * t)
        y = R * np.sin(omega * t)
        theta = np.arctan2(y, x) + np.pi / 2
        kappa = 1.0 / R
        return x, y, theta

    @staticmethod
    def ellipse(t, a=12.0, b=8.0, v=1.0):
        r'''
        Generates an ellipse path and its kinematic properties.
        Returns: $x$, $y$, $\theta$ (heading), $\kappa$ (curvature)
        '''
        omega = v / ((a + b) / 2.0)
        x = a * np.cos(omega * t)
        y = b * np.sin(omega * t)

        dx = -a * omega * np.sin(omega * t)
        dy = b * omega * np.cos(omega * t)
        ddx = -a * omega**2 * np.cos(omega * t)
        ddy = -b * omega**2 * np.sin(omega * t)

        theta = np.arctan2(dy, dx)
        speed_sq = dx**2 + dy**2

        if speed_sq < 1e-6:
            kappa = 0
        else:
            kappa = (dx * ddy - dy * ddx) / speed_sq**(1.5)
        return x, y, theta

    @staticmethod
    def figure_eight(t, a=10.0, v=1.0):
        r'''
        Generates a figure-eight path and its kinematic properties.
        Returns: $x$, $y$, $\theta$ (heading), $\kappa$ (curvature)
        '''
        omega = v / a

        dx = a * omega * np.cos(omega * t)
        dy = a * omega * np.cos(2 * omega * t)
        ddx = -a * omega**2 * np.sin(omega * t)
        ddy = -2 * a * omega**2 * np.sin(2 * omega * t)


        x = a * np.sin(omega * t)
        y = a * np.sin(omega * t) * np.cos(omega * t)
        theta = np.arctan2(dy, dx)
        speed_sq = dx**2 + dy**2

        if speed_sq < 1e-6:
            kappa = 0
        else:
            kappa = (dx * ddy - dy * ddx) / speed_sq**(1.5)
        return x, y, theta

    @staticmethod
    def straight_line(t, slope=0.0, intercept=0.0, v=1.0):
        r'''
        Generates a straight line path and its kinematic properties.
        Returns: $x$, $y$, $\theta$ (heading), $\kappa$ (curvature)
        '''
        x = v * t
        y = slope * x + intercept
        theta = np.arctan2(slope, 1.0)
        kappa = 0.0
        return x, y, theta

    def equation(self, t):
        try:
            if self.shape == "circle":
                return Path.circle(t, R=self.R, v=self.v)
            elif self.shape == "ellipse":
                return Path.ellipse(t, a=self.a, b=self.b, v=self.v)
            elif self.shape == "figure_eight":
                return Path.figure_eight(t, a=self.a, v=self.v)
            elif self.shape == "straight_line":
                return Path.straight_line(t, slope=self.slope, intercept=self.intercept, v=self.v)
        except:
            print("Provide correct path format!")

    def add_obstacle(self, *obstacles):
        '''
        Adds obstacle
        '''
        self.obstacles = obstacles
        pass
