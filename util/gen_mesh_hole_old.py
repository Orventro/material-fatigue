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

def shape_alpha_alter(a,k):
	return 1/(1-a**k/2)-1

def blend(theta, f1, r1, f2, r2, alpha):
	a1 = shape_alpha_alter(alpha, 4)
	return polar_to_decart(theta, (f1(theta)**(1-a1) * f2(theta)**a1) * (r1**(1-alpha) * r2**alpha))

def gen_points(n, m, circle_radius, square_radius):
	ans = []
	theta = np.linspace(0, 2*np.pi, n+1)[:-1] + np.pi/4
	elsizemul = (square_radius**(1/m)*circle_radius**(-1/m) - 1) 
	for alpha in np.linspace(0, 1, m):
		ans.append(blend(theta, lambda x: 1, circle_radius, 
					  	 polar_square, square_radius, alpha))
		if alpha != 0 and alpha != 1:
			elsize = square_radius**alpha*circle_radius**(1-alpha)*elsizemul
			# ans[-1] += np.random.rand(*ans[-1].shape) * elsize * 0.35
			# ans[-1] += np.random.rand(*ans[-1].shape) * elsize * 0.35 * 1j
	return np.concatenate(ans)

def gen_polys(n, m):
	ans = []
	for i in range(m-1):
		f = i*n
		ans.append([f, f+n-1, f+2*n-1, f+n])
		for j in range(0,n-1):
			f = i*n+j
			ans.append([f+1, f, f+n, f+n+1])
	return ans

def gen_mesh(n, m, rc, rs):
	pts = gen_points(n, m, rc, rs)
	polys = gen_polys(n, m)
	boundary = []
	for i in range(n):
		boundary.append((5, 1, (i+1)%n, i))
	for i in range(n):
		boundary.append((int((i)*4/n)%4+1, 1, n*(m-1)+i, n*(m-1)+(i+1)%n))
	return pts, polys, boundary

def print_mesh(pts, polys, n, m, material=None):
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
	if material is None:
		material = [1]*len(polys)
		material = range(1, 1+len(polys))
	print(f'elements\n{len(polys)}')
	for mat, poly in zip(material, polys):
		print(mat, 3, *poly)
	
	print(f'\nboundary\n{n*2}')
	for i in range(n):
		print(5, 1, (i+1)%n, i)
	for i in range(n):
		print(int(i*4/n+3.5)%4+1, 1, n*(m-1)+i, n*(m-1)+(i+1)%n)
	
	print(f'\nvertices\n{len(pts)}\n2')
	for p in pts:
		print(f'{p.real} {p.imag}')

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print('Need 2 arguments: number of points in circle and number of levels')
		exit(1)
	
	n, m = int(sys.argv[1]), int(sys.argv[2])
	
	pts = gen_points(n, m, 0.25e-3, 6e-3)
	
	polys = gen_polys(n, m)

	# material = assign_material(pts, polys, material_whole)
	material = np.arange(len(polys), dtype=int)+1
	
	polys = np.array(polys)

	print_mesh(pts, polys, n, m, material)
	
	exit()
	fig, ax = plt.subplots()
	for p, m in zip(polys, material):
		col = 'b'
		if m == 2:
			col = 'r'
		poly = Polygon(list(zip(np.real(pts[p]), np.imag(pts[p]))), facecolor=col)
		ax.add_patch(poly)
	plt.show()