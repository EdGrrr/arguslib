from cv2 import undistort
import numpy as np


def polynom(x, degree: int):
    """
    Raise `x` with polynomial degrees
    """
    return np.stack(
        [np.ones_like(x), x] + [np.power(x, i) for i in range(2, degree + 1)], axis=-1
    )


def polyval(x, poly):
    """
    Evaluate polynom `poly` at `x`
    """
    return dot(polynom(x, len(poly) - 1), poly)


def image_to_view(poly_rz, principal_point, p_image, norm: bool = False):
    xy = p_image[..., :2] - principal_point
    radius = np.linalg.norm(xy, axis=-1, keepdims=True)
    v_view = np.concatenate((xy, polyval(radius, poly_rz)), axis=-1)
    return unit(v_view) if norm else v_view


def dot(a, b, keepdims=False):
    """
    Dot product in last axis
    """
    return (a * b).sum(axis=-1, keepdims=keepdims)


def unit(data, axis=-1):
    # TODO: Fix for multiple points
    return data / np.linalg.norm(data, axis=axis, keepdims=True)


def get_theta(p, normed: bool = False):
    """
    Angle between points `p` and `(0, 0, 1)`
    """
    if not normed:
        p = unit(p)
    return np.arccos(dot(p, np.array((0, 0, 1)), True))


def project_poly_thetar(
    view_points, poly_theta, principal_point=None, normed: bool = False
):
    """
    Project `view_points` to image points with polynom
    that projects from incident angle`theta` to radius
    `(x, y, z)` is the view vector with angle`theta` to `(0, 0, 1)`
    radius is the length of `(x, y)` vector
    """
    theta = get_theta(view_points, normed)
    rho = polyval(theta, poly_theta)
    xy = unit(view_points[..., :2]) * rho
    if principal_point is not None:
        xy = xy + principal_point
    return xy


class Projection:
    def __init__(self, poly_thetar, poly_rz, principal_point):
        self.poly_thetar = poly_thetar
        self.poly_rz = poly_rz
        self.principal_point = principal_point

    @classmethod
    def fromfile(cls, filename):
        with open(filename) as f:
            import yaml

            calibration = yaml.load(f, Loader=yaml.SafeLoader)

        return cls(
            np.array(calibration["poly_incident_angle_to_radius"]),
            np.array(calibration["poly_radius_to_z"]),
            np.array(calibration["principal_point"]),
        )

    def image_to_view(self, p_image, norm: bool = False):
        """Convert an image pixel location to a 3D location"""
        p_image = np.array(p_image)
        return image_to_view(self.poly_rz, self.principal_point, p_image, norm)

    def view_to_image(self, v_view, normed: bool = False):
        """Convert a 3D location to a pixel location"""
        return project_poly_thetar(
            v_view, self.poly_thetar, self.principal_point, normed
        )


def focal_length_px(focal_length_mm, image_size_px, sensor_size_mm):
    """
    Convert focal length in mm to pixels
    """
    return focal_length_mm * image_size_px / sensor_size_mm


class PerspectiveProjection:
    def __init__(self, focal_lengths, principal_point, distortion_coeffs=None):
        self.focal_lengths = focal_lengths
        self.principal_point = principal_point  # 1x2 containing cx, cy
        self.camera_matrix = np.array(
            [
                [focal_lengths[0], 0, principal_point[0]],
                [0, focal_lengths[1], principal_point[1]],
                [0, 0, 1],
            ]
        )
        self.distortion_coeffs = distortion_coeffs
        # 1x5 containing k1, k2, p1, p2, k3

        if isinstance(distortion_coeffs, list):
            self.distortion_coeffs = np.array(distortion_coeffs)

    def image_to_view(self, p_image, norm: bool = False):
        """
        Convert an image pixel location to a 3D location using OpenCV undistort.
        """
        undistorted = undistort(
            np.expand_dims(p_image, axis=0), self.camera_matrix, self.distortion_coeffs
        )[0]
        xy = undistorted - self.principal_point
        z = np.ones_like(xy[..., :1])  # Assume z = 1 for normalized view direction

        # Standard CV: x=right, y=down, z=forward
        # Arguslib Instrument: x=right, y=up, z=forward

        # CV Frame: X_cv, Y_cv, Z_cv (Forward)
        # Inst Frame: X_inst, Y_inst (Up), Z_inst (Forward)

        # We construct the vector in Instrument frame:
        # X_inst = X_cv
        # Y_inst = -Y_cv (Up in image is -Y in pixels usually)
        # Z_inst = Z_cv (Forward)

        v_view_cv = np.concatenate((xy, z), axis=-1)

        v_view_inst = np.stack(
            [v_view_cv[..., 0], -v_view_cv[..., 1], v_view_cv[..., 2]], axis=-1
        )

        return unit(v_view_inst) if norm else v_view_inst

    def view_to_image(self, v_view, normed: bool = False):
        """
        Convert a 3D location to a pixel location using OpenCV projection.
        Input v_view is in Instrument Frame (X=Right, Y=Up, Z=Forward).
        """
        # Convert Instrument Frame to CV Frame for projection
        # Inst(X, Y, Z) -> CV(X, -Y, Z)
        # X_cv = X_inst
        # Y_cv = -Y_inst
        # Z_cv = Z_inst

        x = v_view[..., 0]
        y = v_view[..., 1]
        z = v_view[..., 2]

        # Perspective division: x' = X_cv / Z_cv, y' = Y_cv / Z_cv
        # x' = x / z
        # y' = -y / z

        # Avoid division by zero
        denom = z
        denom = np.where(np.abs(denom) < 1e-6, 1e-6, denom)

        xy_proj = np.stack([x / denom, -y / denom], axis=-1)

        distorted = xy_proj * self.focal_lengths + self.principal_point
        return distorted
