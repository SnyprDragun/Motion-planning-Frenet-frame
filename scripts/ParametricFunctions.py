#!/Users/subhodeep/venv/bin/python

import numpy as np

class ParametricFunctions:
    """
    Generates parametric paths. Now also returns the path's curvature.
    """
    @staticmethod
    def circle(t, R=10.0, v=1.0):
        """
        Generates a circle path and its kinematic properties.
        Returns: x, y, theta (heading), kappa (curvature)
        """
        omega = v / R
        x = R * np.cos(omega * t)
        y = R * np.sin(omega * t)
        theta = np.arctan2(y, x) + np.pi / 2
        kappa = 1.0 / R
        return x, y, theta, kappa

    @staticmethod
    def ellipse(t, a=12.0, b=8.0, v=1.0):
        """
        Generates an ellipse path and its kinematic properties.
        Returns: x, y, theta (heading), kappa (curvature)
        """
        # Approximate constant velocity
        omega = v / ((a + b) / 2.0)
        
        # Position
        x = a * np.cos(omega * t)
        y = b * np.sin(omega * t)
        
        # First derivatives
        dx = -a * omega * np.sin(omega * t)
        dy = b * omega * np.cos(omega * t)
        
        # Second derivatives
        ddx = -a * omega**2 * np.cos(omega * t)
        ddy = -b * omega**2 * np.sin(omega * t)
        
        # Heading
        theta = np.arctan2(dy, dx)
        
        # Curvature
        speed_sq = dx**2 + dy**2
        if speed_sq < 1e-6:
            kappa = 0
        else:
            kappa = (dx * ddy - dy * ddx) / speed_sq**(1.5)
            
        return x, y, theta, kappa

    @staticmethod
    def figure_eight(t, a=10.0, v=1.0):
        """
        Generates a figure-eight path and its kinematic properties.
        Returns: x, y, theta (heading), kappa (curvature)
        """
        omega = v / a
        
        # First derivatives (velocity components)
        dx = a * omega * np.cos(omega * t)
        dy = a * omega * np.cos(2 * omega * t)
        
        # Second derivatives (acceleration components)
        ddx = -a * omega**2 * np.sin(omega * t)
        ddy = -2 * a * omega**2 * np.sin(2 * omega * t)

        # Path properties
        x = a * np.sin(omega * t)
        y = a * np.sin(omega * t) * np.cos(omega * t)
        theta = np.arctan2(dy, dx)
        
        speed_sq = dx**2 + dy**2
        if speed_sq < 1e-6:
            kappa = 0
        else:
            kappa = (dx * ddy - dy * ddx) / speed_sq**(1.5)

        return x, y, theta, kappa

    @staticmethod
    def straight_line(t, slope=0.0, intercept=0.0, v=1.0):
        """
        Generates a straight line path and its kinematic properties.
        Returns: x, y, theta (heading), kappa (curvature)
        """
        x = v * t
        y = slope * x + intercept
        theta = np.arctan2(slope, 1.0)
        kappa = 0.0
        return x, y, theta, kappa
