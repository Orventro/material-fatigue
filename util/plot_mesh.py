import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import gen_mesh
import sys

def plot_mesh(pts, polys, material=None, filename=None, vmin=None, vmax=None):
    if material is None:
        material = np.zeros(len(polys))
    fig, ax = plt.subplots(1,1, figsize=(10,10))

    if vmin is None:
        vmin = material.min()
    if vmax is None:
        vmax = material.max()

    print(f"vmax={vmax}, vmin={vmin}")

    for poly, mat in zip(polys, material):
        v = (mat - vmin) / (vmax - vmin)
        ax.fill(np.real(pts[poly]), np.imag(pts[poly]), color=[v,v,v])
    plt.axis('off')
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename, bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Need 2 arguments: number of points in circle and number of levels. Material is passed through standard input')
        exit(1)

    n, m = int(sys.argv[1]), int(sys.argv[2])
    pts = gen_mesh.gen_points(n, m, 0.5e-3, 6e-3)
    polys = gen_mesh.gen_polys(n, m)
    sys.stdin.readline() # read psi vector size
    material = np.array(list(map(float, sys.stdin.readlines())))
    # material = np.log(material + 1e-30)

    if len(material) != len(polys):
        print(f'Material input (len={len(material)}) does not conform with provided dimensions ({n}x{m}, {len(polys)} polygons)')
        exit(1)

    plot_mesh(pts, polys, material, None if len(sys.argv) < 4 else sys.argv[3])

    # for j in range(len(polys)):
    #     xs = np.real(pts[polys[j]])
    #     ys = np.imag(pts[polys[j]])
    #     poly = Polygon(np.stack([xs, ys], -1))
    #     # plt.plot(np.real(pts[polys[j]]), np.imag(pts[polys[j]]))
    #     plt.fill(xs, ys, "#afafaf")
    # plt.show()
