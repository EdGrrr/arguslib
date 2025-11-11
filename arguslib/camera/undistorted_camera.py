from .camera import Camera

import numpy as np
import cv2


def _create_remap_arrays(
    output_shape,
    poly_incident_angle_to_radius,
    principal_point,
    focal_length=None,
):
    """
    Pre-computes the x and y remap arrays for undistortion.

    output_shape must be (H, W).
    """
    oh, ow = output_shape  # (height, width)
    cx, cy = map(float, principal_point)

    # Guess focal length if not provided
    if focal_length is None:
        focal_length = abs(float(poly_incident_angle_to_radius[1]))

    # Grid in undistorted space, explicit xy indexing
    x_idx = np.arange(ow, dtype=np.float64)
    y_idx = np.arange(oh, dtype=np.float64)
    i, j = np.meshgrid(x_idx, y_idx, indexing="xy")  # i: x, j: y

    # Normalized rectified plane coords
    x = (i - cx) / focal_length
    y = (j - cy) / focal_length

    # Unit direction vectors
    rays = np.stack([x, y, np.ones_like(x)], axis=-1)
    rays /= np.linalg.norm(rays, axis=-1, keepdims=True)

    # theta and rho
    theta = np.arccos(np.clip(rays[..., 2], -1.0, 1.0))
    rho = np.polyval(poly_incident_angle_to_radius[::-1], theta)

    # Back to distorted pixels
    ray_xy_norm = np.linalg.norm(rays[..., :2], axis=-1, keepdims=True)
    ray_xy_unit = np.divide(
        rays[..., :2],
        ray_xy_norm,
        out=np.zeros_like(rays[..., :2]),
        where=(ray_xy_norm != 0),
    )
    distorted_xy = ray_xy_unit * rho[..., np.newaxis]
    distorted_xy += np.array([cx, cy], dtype=np.float64)

    map_x = distorted_xy[..., 0].astype(np.float32)
    map_y = distorted_xy[..., 1].astype(np.float32)

    # Sanity: maps must be (H, W)
    assert map_x.shape == (oh, ow) and map_y.shape == (oh, ow)
    return map_x, map_y


class UndistortedCamera(Camera):
    """A specialized Camera that applies a fisheye undistortion to images.

    This version pre-computes the undistortion mapping for high performance.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.poly_thetar = self.intrinsic.poly_thetar
        self.intrinsic = UndistortedProjection.from_projection_intrinsics(
            self.intrinsic
        )
        self.reverse_y = False
        self.reverse_x = False

        # Build remap arrays with (H, W)
        W, H = map(int, self.image_size_px)  # self.image_size_px is [W, H]
        output_shape = (H, W)
        self._remap_x, self._remap_y = _create_remap_arrays(
            output_shape,
            self.poly_thetar,
            self.intrinsic.principal_point,
            self.intrinsic.focal_length,
        )

    def get_data_time(self, *args, **kwargs):
        output = super().get_data_time(*args, **kwargs)

        if kwargs.get("return_timestamp", False):
            img, timestamp = output
            undistorted_img = cv2.remap(
                img,
                self._remap_x,
                self._remap_y,
                interpolation=cv2.INTER_LINEAR,
                borderMode=cv2.BORDER_CONSTANT,
            )
            return undistorted_img, timestamp
        else:
            return cv2.remap(
                output,
                self._remap_x,
                self._remap_y,
                interpolation=cv2.INTER_LINEAR,
                borderMode=cv2.BORDER_CONSTANT,
            )

    def target_pix(self, target_position):
        test = super().target_pix(target_position)
        if self.reverse_y:
            test = np.array([test[0], -1 * test[1]])
        if self.reverse_x:
            test = np.array([-1 * test[0], test[1]])
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
