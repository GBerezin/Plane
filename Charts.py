import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

plt.style.use('seaborn-whitegrid')


def geom(coord, joints, npl, nj, name):
    """Геометрия."""

    coord = np.append(coord, coord[0:], axis=0)
    x = coord[:, 0]
    y = coord[:, 1]
    vert = joints[:, :2]
    xj = vert[:, 0]
    yj = vert[:, 1]
    fig, ax = plt.subplots(num=name)
    ax.set_aspect('equal')
    codes = []
    codes += [Path.MOVETO] + [Path.LINETO] * 2 + [Path.CLOSEPOLY]
    ax.scatter(x, y, s=20, c='blue', alpha=0.5)
    for i in range(nj):
        ax.annotate(str(i), (xj[i], yj[i]), size=12, xytext=(
            0, 0), ha='right', c='blue', textcoords='offset points')
    for i in range(npl):
        p3 = coord[i * 3:i * 3 + 3, :]
        vertices = np.append(p3, [[0, 0]], axis=0)
        path = Path(vertices, codes)
        pathpatch = PathPatch(path, facecolor='None', edgecolor='green')
        pcgx = (p3[0, 0] + p3[1, 0] + p3[2, 0]) / 3
        pcgy = (p3[0, 1] + p3[1, 1] + p3[2, 1]) / 3
        ax.add_patch(pathpatch)
        ax.annotate(str(i), (pcgx, pcgy), size=12, xytext=(
            0, 0), ha='center', c='purple', textcoords='offset points')

    plt.title(name, pad=20)
    ax.set_xlabel('X, м')
    ax.set_ylabel('Y, м')
    plt.show()


if __name__ == '__main__':
    print(geom().__doc__)
    input('Press Enter:')
