from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
from .tcrn import model


tspan = np.arange(0, 20*60)
print("Simulating...")
yfull = ScipyOdeSimulator(model).run(tspan=tspan).all

tspan_min = tspan / 60
plt.plot(tspan_min, yfull['A_u'], label='A')
plt.plot(tspan_min , yfull['B_u'], label='B')
plt.plot(tspan_min, yfull['C_u'], label='C')
plt.ylabel("Concentration (uM)")
plt.xlabel("Time (min.)")
plt.legend(loc='upper left')
plt.show()