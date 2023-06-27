"""
Creating mesh for the 2D plate 12 mm x 12 mm with a hole 1 mm in diameter, located in the centre.
"""

import numpy as np
import matplotlib.pylab as plt
from matplotlib.patches import Polygon
import sys
import os

from gen_mesh_hole_old import gen_mesh as gen_mesh_hole_old
from gen_mesh_hole_new import gen_mesh as gen_mesh_hole_new
from gen_mesh_solid import gen_mesh as gen_mesh_solid

def print_mesh(pts, polys, boundary, file):
	file.write("""MFEM mesh v1.0

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
	material = np.arange(len(polys), dtype=int)+1

	file.write(f'elements\n{len(polys)}\n')
	for mat, poly in zip(material, polys):
		file.write(f'{mat} 3 ' + ' '.join(map(str, poly)) + '\n')
	
	file.write(f'\nboundary\n{len(boundary)}\n')
	for i in range(len(boundary)):
		file.write(' '.join(map(str, boundary[i])) + '\n')
	
	file.write(f'\nvertices\n{len(pts)}\n2\n')
	for p in pts:
		file.write(f'{p.real} {p.imag}\n')

if __name__ == '__main__':
	if len(sys.argv) not in [4, 5]:
		print('Need at least 4 arguments: mesh type, number of points in circle, number of levels')
		print('mesh name is optional')
		exit(1)

	mesh_type = sys.argv[1]
	
	n, m = int(sys.argv[2]), int(sys.argv[3])

	rc, rs = 0.5e-3, 1.5e-3

	meshes = {'hole_old' : gen_mesh_hole_old,
			  'hole_new' : gen_mesh_hole_new,
			  'solid' 	 : gen_mesh_solid}


	mesh_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 
												'..', 'meshes'))
	if len(sys.argv) < 5:
		mesh_path = f'{mesh_dir}/mesh_{n}_{m}_{mesh_type}.mesh'
	else:
		mesh_path = f'{mesh_dir}/mesh_{n}_{m}_{sys.argv[4]}.mesh'

	if mesh_type in meshes:
		file = open(mesh_path, 'w')
		pts, polys, boundary = meshes[mesh_type](n, m, rc, rs)
		print_mesh(pts, polys, boundary, file)
		file.close()
		print('Saved to', mesh_path)
	else:
		print(f"Mesh type '{mesh_type}' not recognized.")
		print('Possible values: ' + ', '.join(meshes.keys()))

	