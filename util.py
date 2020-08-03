import openmesh as om
import numpy as np


# center data around origin([0,0,0])
def center_data(mesh):
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

def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix


def rotate_mesh_flat(mesh):
    """ Align directions of greatest variance with x and y axis respectivly
    :param mesh: 3D model
    :return: 3D model rotated to align biggest eigenvector with x axis and second biggest eigenvector with y axis
    """
    point_array = mesh.points()
    cov = np.cov(np.transpose(point_array))
    v, eig = np.linalg.eig(cov)
    print("eigenvalues and vectors")
    print(v)
    print(eig)
    transfMat = rotation_matrix_from_vectors(eig[:, 0],[1,0,0])
    print(transfMat)
    for x in range(len(point_array)):
        point_array[x] = np.matmul(transfMat, point_array[x])

    point_array = mesh.points()
    cov = np.cov(np.transpose(point_array))
    v, eig = np.linalg.eig(cov)
    print("eigenvalues and vectors")
    print(v)
    print(eig)
    transfMat2 = rotation_matrix_from_vectors(eig[:, 1],[0,-1,0])
    for x in range(len(point_array)):
        point_array[x] = np.matmul(transfMat2, point_array[x])
