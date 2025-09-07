#!/Users/subhodeep/venv/bin/python
'''class for conversions'''
import numpy as np

class CartesianFrenetConverter:
    @staticmethod
    def cartesian_to_frenet():
        pass

    @staticmethod
    def frenet_to_cartesian():
        pass

    @ staticmethod
    def normalize_angle(angle):
        '''
        Normalize angle to [-pi, pi]
        '''
        return (angle + np.pi) % (2 * np.pi) - np.pi
