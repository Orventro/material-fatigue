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
	theta = np.linspace(0, 2*np.pi, n+1)[:-1]
	for alpha in np.linspace(0, 1, m):
		ans.append(blend(theta, lambda x: 1, circle_radius, 
					  	 polar_square, square_radius, alpha))
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

def material_eye(pt):
	x = np.real(pt) * 400
	y = np.imag(pt) * 400
	return ((x**2 + y**2) * np.exp(2*np.abs(y)**.5 - 1) <= 1) + 1

def material_circle(pt):
	x = np.real(pt) * 1e3
	y = np.imag(pt) * 1e3
	return ((x**2 + y**2) <= 4) + 1

def material_whole(pt):
	return 1

def assign_material(pts, polys, formula):
	mid = np.mean(pts[np.array(polys)], axis=1) # middle of polygon
	return formula(mid)

def print_vtk(pts, polys, material=None):
	print('# vtk DataFile Version 3.0')
	print('Generated by a script')
	print('ASCII')
	print('DATASET UNSTRUCTURED_GRID')

	print(f'POINTS {len(pts)} double')
	for p in pts:
		print(f'{p.real} {p.imag} 0')

	print(f'CELLS {len(polys)} {5*len(polys)}')
	for poly in polys:
		print(4, *poly)
	
	print(f'CELL_TYPES {len(polys)}')
	print('\n'.join(['9']*len(polys)))

	print(f'CELL_DATA {len(polys)}')
	print(f'SCALARS material int')
	print(f'LOOKUP_TABLE default')
	if material is None:
		print('\n'.join(map(str, [1]*len(polys))))
	else:
		print('\n'.join(map(str, material)))

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
	
	pts = gen_points(n, m, 0.5e-3, 6e-3)
	
	polys = gen_polys(n, m)

	material = assign_material(pts, polys, material_circle)
	
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