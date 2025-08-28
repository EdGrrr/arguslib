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

    This function contains the expensive, one-off calculations.
    """
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
    # Handle the case where rays are straight ahead to avoid division by zero
    ray_xy_norm = np.linalg.norm(rays[..., :2], axis=-1, keepdims=True)
    ray_xy_unit = np.divide(rays[..., :2], ray_xy_norm, out=np.zeros_like(rays[...,:2]), where=(ray_xy_norm != 0))

    distorted_xy = ray_xy_unit * rho[..., np.newaxis]
    distorted_xy += np.array(principal_point)

    # Return the remap arrays
    map_x = distorted_xy[..., 0].astype(np.float32)
    map_y = distorted_xy[..., 1].astype(np.float32)

    return map_x, map_y

class UndistortedCamera(Camera):
    """A specialized Camera that applies a fisheye undistortion to images.

    This version pre-computes the undistortion mapping for high performance.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.poly_thetar = self.intrinsic.poly_thetar
        self.intrinsic = UndistortedProjection.from_projection_intrinsics(self.intrinsic)
        self.reverse_y=False
        self.reverse_x=False

        # --- OPTIMIZATION ---
        # Pre-compute the remap arrays. We assume the output shape is the same
        # as the camera's native shape. You might need to adjust self.shape.
        # For example, use self.intrinsic.image_shape if available.
        output_shape = self.image_size_px # Assumes self.shape exists from the base Camera class
        focal_length = self.intrinsic.focal_length
        principal_point = self.intrinsic.principal_point
        
        self._remap_x, self._remap_y = _create_remap_arrays(
            output_shape, self.poly_thetar, principal_point, focal_length
        )
    
    def get_data_time(self, *args, **kwargs):
        output = super().get_data_time(*args, **kwargs)
        
        if kwargs.get('return_timestamp', False):
            img, timestamp = output
            undistorted_img = cv2.remap(img, self._remap_x, self._remap_y, interpolation=cv2.INTER_LINEAR, borderMode=cv2.BORDER_CONSTANT)
            return undistorted_img, timestamp
        else:
            return cv2.remap(output, self._remap_x, self._remap_y, interpolation=cv2.INTER_LINEAR, borderMode=cv2.BORDER_CONSTANT)
    
    def target_pix(self, target_position):
        test = super().target_pix(target_position)
        if self.reverse_y:
            test = np.array([test[0], -1* test[1]])
        if self.reverse_x:
            test = np.array([-1* test[0], test[1]])
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