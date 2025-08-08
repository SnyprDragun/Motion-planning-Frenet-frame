#!/Users/subhodeep/venv/bin/python

import numpy as np

class ParametricFunctions:
    @staticmethod
    def circle(t, R=10.0, v=1.0):
        omega = v / R
        x = R * np.cos(omega * t)
        y = R * np.sin(omega * t)
        theta = np.arctan2(y, x) + np.pi / 2
        return x, y, theta

    @staticmethod
    def ellipse(t, a=12.0, b=8.0, v=1.0):
        omega = v / ((a + b) / 2.0)
        x = a * np.cos(omega * t)
        y = b * np.sin(omega * t)
        theta = np.arctan2(b * np.cos(omega * t), -a * np.sin(omega * t))
        return x, y, theta

    @staticmethod
    def figure_eight(t, a=10.0, v=1.0):
        omega = v / a
        x = a * np.sin(omega * t)
        y = a * np.sin(omega * t) * np.cos(omega * t)
        dx = a * omega * np.cos(omega * t)
        dy = a * omega * np.cos(2 * omega * t)
        theta = np.arctan2(dy, dx)
        return x, y, theta
