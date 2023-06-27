"""
Creating mesh for the 2D plate 12 mm x 12 mm with a hole 1 mm in diameter, located in the centre.
"""

import numpy as np
import matplotlib.pylab as plt
from matplotlib.patches import Polygon
import sys

def polar_square(theta):
	return np.minimum(1/np.abs(np.sin(theta)+1e-12), 1/np.abs(np.cos(theta)+1e-12))

def polar_to_decart(theta, radius):
	return np.exp(1j*theta) * radius

def shape_alpha_alter(a, k):
	return 1/(1-a**k/2)-1

def blend(theta, f1, r1, f2, r2, alpha):
	a1 = shape_alpha_alter(alpha, 2)
	# return polar_to_decart(theta, (f1(theta)**(1-a1) * f2(theta)**a1) * (r1**(1-alpha) * r2**alpha))
	return polar_to_decart(theta, f1(theta)*(1-alpha**2) + f2(theta)*alpha**2) * (r1 + (r2-r1)*alpha)

def gen_mesh(n=96, m=50, rc=0.5e-3, rs=6e-3):
	"""Generate mesh. 
		Arguments:
			n: number of points in the circle
			m: number of layers
			rc: circle radius
			rs: square radius
		Returns:
			points, polygons, polygon types, boundary
	"""
	pts = []
	polys = []
	boundary = []
	f = 0
	N = n
	theta = np.linspace(0, 2*np.pi, N+1)[:-1] + np.pi/4
	layer = polar_to_decart(theta, rc)
	for i in range(N):
		boundary.append([5, 1, i, (i+1)%N])
	pts.append(layer)
	for i in range(1,m):
		gamma = (rc*(1-i/(m-1)) + rs*i/(m-1))/rc
		if 4*N <= gamma*n:
			polys.append([f, f+N*3-1, f+N, f+N+1])
			polys.append([f, f+N+1, f+N+2, f+1])
			for j in range(2, N, 2):
				polys.append([f+j, f+j-1, f+j*2+N-2, f+j*2+N-1])
				polys.append([f+j, f+j*2+N-1, f+j*2+N, f+j*2+N+1])
				polys.append([f+j, f+j*2+N+1, f+j*2+N+2, f+j+1])
			polys.append([f, f+N-1, f+N*3-2, f+N*3-1])
			f += N
			N *= 2
			theta = np.linspace(0, 2*np.pi, N+1)[:-1] + np.pi/4
			layer = blend(theta, lambda x: 1, rc, polar_square, rs, i/(m-1))
			layer *= 1+(np.arange(1, N+1, dtype=int)%2*2-1)*0.25/m
		else:
			for j in range(0,N-1):
				polys.append([f+j, f+j+N, f+j+N+1, f+j+1])
			polys.append([f, f+N-1, f+2*N-1, f+N])
			f += N
			layer = blend(theta, lambda x: 1, rc, polar_square, rs, i/(m-1))
		pts.append(layer)
	for i in range(N):
		boundary.append([i//(N//4)+1, 1, f+i, f+(i+1)%N])
	return np.concatenate(pts), polys,  np.array(boundary)


def print_mesh(pts, polys, boundary, material):
	print("""MFEM mesh v1.0

#
# MFEM Geometry Types (see mesh/geom.hpp):
#
# POINT       = 0
# SEGMENT     = 1
# TRIANGLE    = 2
# SQUARE      = 3
# TETRAHEDRON = 4
# CUBE        = 5
# PRISM       = 6
#

dimension
2
""")
	print(f'elements\n{len(polys)}')
	for mat, poly in zip(material, polys):
		print(mat, 3, *poly)
	
	print(f'\nboundary\n{boundary.shape[0]}')
	for i in range(boundary.shape[0]):
		print(*boundary[i])
	
	print(f'\nvertices\n{len(pts)}\n2')
	for p in pts:
		print(p.real, p.imag)

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print('Need 2 arguments: number of points in circle and number of levels')
		exit(1)
	
	n, m = int(sys.argv[1]), int(sys.argv[2])

	pts, polys, boundary = gen_mesh(n, m, 0.5e-3, 6e-3)

	# material = assign_material(pts, polys, material_whole)
	material = np.arange(len(polys), dtype=int)+1
	print_mesh(pts, polys, boundary, material)
	
	exit()
	fig, ax = plt.subplots()
	for p, m in zip(polys, material):
		col = 'b'
		if m == 2:
			col = 'r'
		poly = Polygon(list(zip(np.real(pts[p]), np.imag(pts[p]))), facecolor=col)
		ax.add_patch(poly)
	plt.show()