from sympy import Plane, Point, Point3D, Segment3D, Line3D, Ray3D
from sympy.core import Dummy, Rational, S, Symbol
import numpy as np
from numpy.linalg import norm





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
    """
    Align directions of greatest variance with x and y axis respectivly
    Rotate 3D model to align biggest eigenvector with x axis and second biggest eigenvector with y axis
    :param mesh: 3D model
    """
    point_array = mesh.points()
    cov = np.cov(np.transpose(point_array))
    v, eig = np.linalg.eig(cov)
    print("eigenvalues and vectors")
    print(v)
    print(eig)
    transfMat = rotation_matrix_from_vectors(eig[:, 0], [1, 0, 0])
    print(transfMat)
    for x in range(len(point_array)):
        point_array[x] = np.matmul(transfMat, point_array[x])

    point_array = mesh.points()
    cov = np.cov(np.transpose(point_array))
    v, eig = np.linalg.eig(cov)
    print("eigenvalues and vectors")
    print(v)
    print(eig)
    transfMat2 = rotation_matrix_from_vectors(eig[:, 1], [0, -1, 0])
    for x in range(len(point_array)):
        point_array[x] = np.matmul(transfMat2, point_array[x])


def edge_is_same_side(plane, edge):
    """ Returns true if both points of given edge are on the same side of the given plane
    :param plane: sympy plane
    :param edge: Two points in form of: [point0,point1] => [(1,1,1),(2,1,2)]
    :return: true or false
    """
    point0 = Point(edge[0])
    point1 = Point(edge[1])
    sign0 = check_side_of_plane(plane, point0)
    sign1 = check_side_of_plane(plane, point1)
    print(sign0)
    print(sign1)
    return sign0 == sign1


def check_side_of_plane(plane, point):
    """ Checks if point is located above, on or under given plane
    :param plane:
    :param point:
    :return: 1 for above plane,0 for on plane or -1 for under plane
    """
    result = distance(plane, point)
    if (result == 0):
        return 0
    sign = abs(result) / result
    return sign


def distance(self, o):
    """ From https://github.com/sympy/sympy/blob/58e1e9abade7c38c66871fd06bf76ebac8c1df78/sympy/geometry/plane.py#L246-L298
        but returns signed Distance not absolut
    :param self: Plane
    :param o: Point
    :return: signed distance between Plane and Point
    """
    if self.intersection(o) != []:
        return S.Zero

    if isinstance(o, (Segment3D, Ray3D)):
        a, b = o.p1, o.p2
        pi, = self.intersection(Line3D(a, b))
        if pi in o:
            return self.distance(pi)
        elif a in Segment3D(pi, b):
            return self.distance(a)
        else:
            assert isinstance(o, Segment3D) is True
            return self.distance(b)

    # following code handles `Point3D`, `LinearEntity3D`, `Plane`
    a = o if isinstance(o, Point3D) else o.p1
    n = Point3D(self.normal_vector).unit
    d = (a - self.p1).dot(n)
    return d
