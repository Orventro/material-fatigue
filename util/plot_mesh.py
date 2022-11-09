import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def plot_mesh(pts, polys, material=None, filename=None, vmin=None, vmax=None):
    if material is None:
        material = np.zeros(len(polys))
    fig, ax = plt.subplots(1,1, figsize=(10,10))

    if vmin is None:
        vmin = material.min()
    if vmax is None:
        vmax = material.max()

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
    input()
    input()
    input()
    input()

    _, pts_num, _ = input().split()
    pts_num = int(pts_num)

    pts = []
    for i in range(pts_num):
        x, y, z = input().split()
        pts.append(float(x)+1j*float(y))
    pts = np.array(pts)

    _, poly_num, _ = input().split()
    poly_num = int(poly_num)

    for j in range(poly_num):
        arr = list(map(int, input().split()))
        arr[0] = arr[-1]
        plt.plot(np.real(pts[arr]), np.imag(pts[arr]))
    plt.show()
