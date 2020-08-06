import numpy as np
import openmesh as om

# center data around origin([0,0,0])
def center_data(mesh):
    """
    Centers mesh data
    :param mesh: mesh to center
    """
    xValues = 0;
    yValues = 0;
    zValues = 0;
    point_array = mesh.points()
    for node in point_array:
        xValues += node[0]
        yValues += node[1]
        zValues += node[2]
    xMean = xValues / len(point_array)
    yMean = yValues / len(point_array)
    zMean = zValues / len(point_array)
    print(len(point_array))
    print(xMean)
    print(yMean)
    print(zMean)
    point_array -= np.array([xMean, yMean, zMean])