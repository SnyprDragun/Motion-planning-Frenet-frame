#!/Users/subhodeep/venv/bin/python

import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("curve_frenet_output.csv")

fig, ax = plt.subplots(1,2,figsize=(12,5))
ax[0].plot(data.x, data.y, 'b')
ax[0].set_title("Curve in Cartesian Frame")
ax[0].set_xlabel("x"); ax[0].set_ylabel("y"); ax[0].axis('equal')

ax[1].plot(data.s, data.d, 'r')
ax[1].set_title("Curve in Frenet Frame")
ax[1].set_xlabel("s (arc length)"); ax[1].set_ylabel("d (lateral offset)")

plt.tight_layout(); plt.show()
