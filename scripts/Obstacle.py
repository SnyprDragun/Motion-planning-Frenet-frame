#!/Users/subhodeep/venv/bin/python

from dataclasses import dataclass

@dataclass
class Obstacle:
    x: float
    y: float
    radius: float
