import numpy as np

with open('cam1_calibration.yml') as f:
    import yaml
    calibration = yaml.load(f, Loader=yaml.SafeLoader)

poly_rz = np.array(calibration['poly_radius_to_z'], dtype=np.float32)
poly_thetar = np.array(calibration['poly_incident_angle_to_radius'])
principal_point = np.array(calibration['principal_point'])

def polynom(x, degree: int):
    '''
    Raise `x` with polynomial degrees
    '''
    return np.stack([
        np.ones_like(x), x] + [np.power(x, i)
                               for i in range(2, degree + 1)
                               ], axis=-1)


def polyval(x, poly):
    '''
    Evaluate polynom `poly` at `x`
    '''
    return dot(polynom(x, len(poly) - 1), poly)


def image_to_view(poly_rz, principal_point, p_image,
                  norm: bool = False):
    xy = p_image[..., :2] - principal_point
    radius = np.linalg.norm(
        xy,
        axis=-1,
        keepdims=True
    )
    v_view = np.concatenate((
        xy,
        polyval(radius, poly_rz)
    ), axis=-1)
    return unit(v_view) if norm else v_view


def dot(a, b, keepdims=False):
    '''
    Dot product in last axis
    '''
    return (a * b).sum(axis=-1, keepdims=keepdims)

def unit(data, axis=-1):
    # TODO: Fix for multiple points
    return data/np.linalg.norm(data, axis=axis, keepdims=True)


def get_theta(p, normed: bool = False):
    '''
    Angle between points `p` and `(0, 0, 1)`
    '''
    if not normed:
        p = unit(p)
    return np.arccos(dot(p, np.array((0, 0, 1)), True))

def project_poly_thetar(view_points, poly_theta,
                        principal_point = None,
                        normed: bool = False):
    '''
    Project `view_points` to image points with polynom
    that projects from incident angle`theta` to radius
    `(x, y, z)` is the view vector with angle`theta` to `(0, 0, 1)`
    radius is the length of `(x, y)` vector
    '''
    theta = get_theta(view_points, normed)
    rho = polyval(theta, poly_theta)
    xy = unit(view_points[..., :2]) * rho
    if principal_point is not None:
        xy = xy + principal_point
    return xy


class Projection:
    def __init__(self, poly_thetar, poly_rz,
                 principal_point):
        self.poly_thetar = poly_thetar
        self.poly_rz = poly_rz
        self.principal_point = principal_point

    @classmethod
    def fromfile(cls, filename):
        with open(filename) as f:
            import yaml
            calibration = yaml.load(f, Loader=yaml.SafeLoader)

        return cls(np.array(calibration['poly_incident_angle_to_radius']),
                   np.array(calibration['poly_radius_to_z']),
                   np.array(calibration['principal_point']))

    def image_to_view(self, p_image, norm: bool = False):
        '''Convert an image pixel location to a 3D location'''
        return image_to_view(
            self.poly_rz,
            self.principal_point,
            p_image,
            norm
        )

    def view_to_image(self, v_view, normed: bool = False):
        '''Convert a 3D location to a pixel location'''
        return project_poly_thetar(
            v_view,
            self.poly_thetar,
            self.principal_point,
            normed
        )

# projection = Projection(poly_thetar, poly_rz, principal_point)
# landmark_p = np.array((2034, 2788))
# landmark_v = projection.image_to_view(landmark_p, True)
# print(landmark_v)
# landmark_v = projection.image_to_view(np.array([landmark_p, landmark_p]), True)
# print(landmark_v)

# print(landmark_p)
# landmark_pn = projection.view_to_image(landmark_v, True)
# print(landmark_pn)

# p2 = Projection.fromfile('cam1_calibration.yml')
# landmark_v = p2.image_to_view(landmark_p, True)
# print(landmark_v)
