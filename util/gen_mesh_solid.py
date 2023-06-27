"""
Creating mesh for the 2D plate 12 mm x 12 mm with a hole 1 mm in diameter, located in the centre.
"""

import numpy as np
import matplotlib.pylab as plt
from matplotlib.patches import Polygon
import sys

def gen_mesh(n, m, rc, rs):
	N = n//4+1
	pts = np.linspace(-rs, rs, N)[None,:] + np.linspace(-rs, rs, N)[:,None]*1j
	pts = pts.flatten()
	polys = []
	for i in range(N-1):
		for j in range(N-1):
			polys.append((i*N+j, (i+1)*N+j, (i+1)*N+j+1, i*N+j+1))
	boundary = []
	for i in range(N-1):
		boundary.append((1, 1, i+1, i))
	for i in range(N-1):
		boundary.append((2, 1, i*N, (i+1)*N))
	for i in range(N-1):
		boundary.append((3, 1, (N-1)*N+i, (N-1)*N+i+1))
	for i in range(N-1):
		boundary.append((4, 1, (i+1)*N+N-1, i*N+N-1))
	return pts, polys, boundary