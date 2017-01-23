import numpy as np
import matplotlib.pyplot as plt 
import os
from scipy.signal import argrelextrema

data = np.load('./output-2016_10_04_23_22_25/overlaps_eta_0.08_chi1_0.50_N_100.npz')

slices = np.arange(0,100,1)
points = []

for slice in slices:
	somedata = data['OLVP_0F_P2'][:,slice]
	index = argrelextrema(somedata, np.greater)[0]
	if len(index) == 1 and index > 60 and index < 80:
		points.append(np.array([slice, index + 1]))

points = np.array(points)
coeffs = np.polyfit(points[:,0], points[:,1], 1)

x = points[:,0]
y = points[:,1]
z = coeffs[0]*x + coeffs[1]

plt.plot(x, y, 'm.')
plt.plot(x, z, 'w-', linewidth=2.0)
plt.imshow(data['OLVP_0F_P2'], cmap='gray_r')
plt.show()