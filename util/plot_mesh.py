import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import gen_mesh
import sys

def plot_mesh(pts, polys, material=None, filename=None, vmin=None, vmax=None):
    if material is None:
        material = np.zeros(len(polys))
    fig, ax = plt.subplots(1,1, figsize=(10,10))

    if vmin is None:
        vmin = material.min()
    else:
        vmin = min(material.min(), vmin)
    if vmax is None:
        vmax = material.max()
    else:
        vmax = max(material.max(), vmax)

    print(f"vmax={vmax}({material.max()}), vmin={vmin}({material.min()})")

    for poly, mat in zip(polys, material):
        v = (mat - vmin) / (vmax - vmin + 1e-12)
        ax.fill(np.real(pts[poly]), np.imag(pts[poly]), color=[v,v,v])
    plt.axis('off')
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename, bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                prog='PlotMesh',
                description='Plot mesh with provided material')
    parser.add_argument('material', help='Material filename')
    parser.add_argument('circ', help='Number of points in mesh base circle', type=int)
    parser.add_argument('nlev', help='Number of levels in mesh', type=int)
    parser.add_argument('-n', '--name', help='Name of column of material')
    parser.add_argument('--min', help='Minimum value of material', type=float)
    parser.add_argument('--max', help='Maximum value of material', type=float)
    parser.add_argument('-o', '--output', help='Output filename')
    args = parser.parse_args()

    print(args)

    pts = gen_mesh.gen_points(args.circ, args.nlev, 0.25e-3, 6e-3)
    polys = gen_mesh.gen_polys(args.circ, args.nlev)

    try:
        material = pd.read_csv(args.material)
        if not args.name:
            material = material.iloc[:,0].to_numpy()
        else:
            try:
                material = material[args.name]
            except Exception as e:
                print('Could not find column', e)
                exit(1)
    except Exception as e:
        print('Could not read material. Here\'s the exception:')
        print(e)
        exit(1)

    if len(material) != len(polys):
        print(f'Material input (len={len(material)}) does not conform with provided dimensions ({n}x{m}, {len(polys)} polygons)')
        exit(1)

    plot_mesh(pts, polys, material, args.output, args.min, args.max)
