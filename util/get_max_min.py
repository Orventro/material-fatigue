import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import gen_mesh
import sys

def get_max_min(fnames, material, mod, silent=True):
    if mod == 'min':
        ans = np.infty
    else:
        ans = -np.infty
    for fn in fnames:
        try:
            df = pd.read_csv(fn)
        except Exception as e:
            if not silent:
                print(f'Could not read file {fn}:', e)
            continue
        try:
            arr = df[material].to_numpy()
        except Exception as e:
            if not silent:
                print(f'File {fn} has no {material} column.', e)
            continue
        if mod == 'min':
            ans = np.minimum(ans, np.nanmin(arr))
        else:
            ans = np.maximum(ans, np.nanmax(arr))
    return ans

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                prog='PlotMesh',
                description='Plot mesh with provided material')
    parser.add_argument('filenames', help='Checked files', nargs="+")
    parser.add_argument('material', help='Material name', nargs=1)
    parser.add_argument('mode', choices=['min', 'max'], help='Return min or max')
    args = parser.parse_args()
    
    print(get_max_min(args.filenames, args.material, args.mode))
