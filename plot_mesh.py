import numpy as np
import matplotlib.pyplot as plt

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
