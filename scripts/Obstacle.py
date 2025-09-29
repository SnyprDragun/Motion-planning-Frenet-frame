#!/Users/subhodeep/venv/bin/python3
'''class to add obstacle'''
from dataclasses import dataclass

@dataclass
class Obstacle:
    x: float
    y: float
    radius: float
