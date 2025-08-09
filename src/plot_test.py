#!/Users/subhodeep/venv/bin/python

import pandas as pd
import matplotlib.pyplot as plt

s = pd.read_csv('states.csv')
h = pd.read_csv('hitches.csv')
t = pd.read_csv('trailers.csv')
f = pd.read_csv('frenet.csv')

plt.figure(figsize=(8,8))
plt.plot(s.x, s.y, '-r', label='Mule Path')
plt.plot(h.x, h.y, 'yo', markersize=3, label='Hitch')
plt.plot(t.x, t.y, 'k', label='Trailer')
plt.legend(); plt.axis('equal'); plt.show()

plt.figure()
plt.plot(f.s, f.d, '-r'); plt.axhline(0, color='gray', linestyle='--')
plt.xlabel('s'); plt.ylabel('d'); plt.show()
