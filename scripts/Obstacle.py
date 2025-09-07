#!/Users/subhodeep/venv/bin/python
'''class to add obstacle'''
from dataclasses import dataclass

@dataclass
class Obstacle:
    x: float
    y: float
    radius: float
