
from .camera import Camera

import numpy as np
import cv2

def undistort_custom_fisheye(
    image,
    poly_incident_angle_to_radius,
    principal_point,
    focal_length=None,
    output_shape=None,
):
    """
    Undistort a fisheye image using a custom θ → ρ polynomial model.

    Parameters:
        image: distorted input image (H, W, 3) or (H, W)
        poly_incident_angle_to_radius: list of polynomial coefficients
        principal_point: (cx, cy) in pixels
        focal_length: optional scalar (fx = fy), default: inferred from poly
        output_shape: (height, width) of undistorted image, default = input image size
    Returns:
        undistorted_image
    """
    h, w = image.shape[:2]
    if output_shape is None:
        output_shape = (h, w)
    oh, ow = output_shape
    cx, cy = principal_point

    # Guess focal length if not provided
    if focal_length is None:
        # Use the linear coefficient as a rough estimate
        focal_length = abs(poly_incident_angle_to_radius[1])

    # Build meshgrid in undistorted (rectified) space
    i, j = np.meshgrid(np.arange(ow), np.arange(oh))
    x = (i - cx) / focal_length
    y = (j - cy) / focal_length

    # Unit direction vectors from optical center
    rays = np.stack([x, y, np.ones_like(x)], axis=-1)
    rays /= np.linalg.norm(rays, axis=-1, keepdims=True)

    # Compute θ = angle from optical axis
    theta = np.arccos(rays[..., 2])

    # Compute ρ = distorted radius via polynomial
    rho = np.polyval(poly_incident_angle_to_radius[::-1], theta)

    # Project back to distorted image coordinates
    ray_xy_unit = rays[..., :2] / np.linalg.norm(rays[..., :2], axis=-1, keepdims=True)
    distorted_xy = ray_xy_unit * rho[..., np.newaxis]
    distorted_xy += np.array(principal_point)

    # Build remap arrays
    map_x = distorted_xy[..., 0].astype(np.float32)
    map_y = distorted_xy[..., 1].astype(np.float32)

    # Remap
    undistorted = cv2.remap(image, map_x, map_y, interpolation=cv2.INTER_LINEAR, borderMode=cv2.BORDER_CONSTANT)

    return undistorted

class UndistortedCamera(Camera):
    """A specialized Camera that applies a fisheye undistortion to images.

    This class inherits from `Camera` and overrides the `get_data_time` method.
    When data is fetched, it applies a custom polynomial undistortion model
    to rectify the fisheye image into a perspective-like view before it is
    returned or displayed. All other `Camera` functionality, such as geolocation
    and annotation, remains the same.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.poly_thetar = self.intrinsic.poly_thetar
        self.intrinsic = UndistortedProjection.from_projection_intrinsics(self.intrinsic)
        self.reverse_y=False
    
    def get_data_time(self, *args, **kwargs):
        output = super().get_data_time(*args, **kwargs)
        if kwargs.get('return_timestamp', False):
            img, timestamp = output
            return undistort_custom_fisheye(img, self.poly_thetar, self.intrinsic.principal_point), timestamp
        else:
            return undistort_custom_fisheye(output, self.poly_thetar, self.intrinsic.principal_point)
    
    def target_pix(self, target_position):
        test = super().target_pix(target_position)
        if self.reverse_y:
            return np.array([test[0], -1* test[1]])
        else:
            return test
        

class UndistortedProjection:
    def __init__(self, focal_length, principal_point):
        self.focal_length = focal_length
        self.principal_point = np.array(principal_point)

    def view_to_image(self, v_view, normed=False):
        x = v_view[..., 0] / v_view[..., 2]
        y = v_view[..., 1] / v_view[..., 2]
        u = self.focal_length * x + self.principal_point[0]
        v = self.focal_length * y + self.principal_point[1]
        return np.stack([u, v], axis=-1)

    def image_to_view(self, p_image, norm: bool = False):
        p_image = np.array(p_image)
        x = (p_image[..., 0] - self.principal_point[0]) / self.focal_length
        y = (p_image[..., 1] - self.principal_point[1]) / self.focal_length
        z = np.ones_like(x)
        vec = np.stack([x, y, z], axis=-1)
        return vec / np.linalg.norm(vec, axis=-1, keepdims=True) if norm else vec
    
    @classmethod
    def from_projection_intrinsics(cls, projection):
        # guess focal length
        # Use the linear coefficient as a rough estimate
        focal_length = abs(projection.poly_thetar[1])
        return cls(focal_length, projection.principal_point)