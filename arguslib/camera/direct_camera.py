from .undistorted_camera import UndistortedCamera
from .camera import Camera
import numpy as np
from typing_extensions import override
from ..misc.geo import haversine
from ..misc.plotting import plot_range_rings


class DirectCamera(
    Camera
):  # FIXME: annotation (range rings) is broken for direct camera but not direct undistorted
    """A Camera subclass that renders annotations directly onto the image array.

    This class is optimized for performance, especially for creating videos.
    Instead of using Matplotlib for plotting, it overrides the `show` and
    `annotate_positions` methods to use OpenCV functions (`cv2.line`, `cv2.putText`)
    to draw directly on the image data.

    As a result, its `show` and `annotate_positions` methods do not accept or
    return Matplotlib `Axes` objects. The final image can be accessed via the
    `.image` property or `.to_image_array()` method.
    """

    def get_data_time(self, *args, **kwargs):
        ret_val = super().get_data_time(*args, **kwargs)
        self.data_loader.image = self.data_loader.image[:, :, ::-1]
        if isinstance(ret_val, tuple):
            return self.data_loader.image, *ret_val[1:]
        else:
            return self.data_loader.image

    @override
    def annotate_positions(
        self, positions, dt, ax, *args, plotting_method=None, max_range_km=90, **kwargs
    ):
        from matplotlib.colors import to_rgb
        import cv2

        if ax is not None:
            raise ValueError("Direct instruments should have no axes!")

        # TODO: this should take dt into account, mostly because the calibration may change for the same camera at different times...
        lats = [p.lat for p in positions]
        lons = [p.lon for p in positions]

        dists = np.array(
            [
                haversine(self.position.lon, self.position.lat, lon, lat)
                for lon, lat in zip(lons, lats)
            ]
        )

        behind_camera = np.array([self.target_iead(p)[0] < 0 for p in positions])

        if (dists[~np.isnan(dists)] > max_range_km).all():
            return

        pl_track = np.array([self.target_pix(p) for p in positions])

        color = (255, 0, 0)
        if (kw_col := kwargs.get("color")) is not None:
            color = tuple((np.array(to_rgb(kw_col)) * 255).astype(int).tolist())
        elif (kw_col := kwargs.get("c")) is not None:
            color = tuple((np.array(to_rgb(kw_col)) * 255).astype(int).tolist())

        linewidth = 3
        if (kw_linewidth := kwargs.get("kw_linewidth")) is not None:
            linewidth *= kw_linewidth
        if (kw_linewidth := kwargs.get("linewidth")) is not None:
            linewidth *= kw_linewidth

        if plotting_method is None:
            for i in range(pl_track.shape[0] - 1):
                if behind_camera[i] or behind_camera[i + 1]:
                    continue
                if np.isnan(pl_track[i : i + 2]).any():
                    continue
                xy0 = tuple(pl_track[i][0:2].astype(int).tolist())
                xy1 = tuple(pl_track[i + 1][0:2].astype(int).tolist())
                # print(xy0, xy1)
                cv2.line(
                    self.data_loader.image, xy0, xy1, color, thickness=int(linewidth)
                )
        else:
            raise ValueError("Online line plotting for DirectCameras")

        return ax

    @override
    def _show(
        self,
        dt,
        ax=None,
        fail_if_no_data=True,
        imshow_kw={},
        brightness_adjust=1.0,
        **kwargs,
    ):
        """Show the nearest possible timestamp.

        Limited by camera time resolution (5s)"""

        if ax is not None:
            raise ValueError("No axes for DirectCameras")
        if imshow_kw != {}:
            raise ValueError("No imshow for DirectCameras")

        defaults = {"theta_behaviour": "bearing", "lr_flip": True}

        if "theta_behaviour" in kwargs and "lr_flip" not in kwargs:
            # if theta_behaviour is set, assume lr_flip should be false, unless it's bearing
            if kwargs["theta_behaviour"] != "bearing":
                defaults["lr_flip"] = False

        if "lr_flip" in kwargs and "theta_behaviour" not in kwargs:
            # if lr_flip is set, but we are inferring theta_behaviour, assume it's bearing if flipped, ordinal aligned if not flipped
            if kwargs["lr_flip"]:
                defaults["theta_behaviour"] = "bearing"
            else:
                defaults["theta_behaviour"] = "unflipped_ordinal_aligned"

        kwargs = defaults | kwargs

        self.direct_lr_flip = kwargs.pop("lr_flip")
        self.direct_theta_behaviour = kwargs.pop("theta_behaviour")

        try:
            img, timestamp = self.get_data_time(dt, return_timestamp=True)
            # ax.get_figure().timestamp = timestamp
        except FileNotFoundError as e:
            if fail_if_no_data:
                raise e
            else:
                return ax

        # new_vmin = np.round(brightness_adjust * 255)
        self.data_loader.image = np.clip(
            np.uint16(img) * brightness_adjust, 0, 255
        ).astype(np.uint8)
        # img[img <= new_vmin] = new_vmin
        # img = np.round(255 * (img.astype(float) - new_vmin) / (255 - new_vmin)).astype(
        #     int
        # )

        plot_range_rings(self, self, dt, ax=ax)
        # ax.imshow(img[:, :, ::-1], origin="upper", **imshow_kw)

        return None  # No ax for image

    @property
    def image(self):
        from PIL import Image

        return Image.fromarray(
            self.to_image_array()
        )  # Watch out! this image looks like the mp4 files. vertical flip to get it as matplotlib might plot.

    def to_image_array(self, time=True):
        import cv2

        slicer = (
            slice(None, None, -1),
            slice(None, None, None if self.direct_lr_flip else -1),
            slice(None, None, None),
        )
        im = self.data_loader.image
        if time:
            cv2.putText(
                im,
                self.data_loader.current_image_time.isoformat(timespec="seconds"),
                (50, 2900),
                1,
                10,
                (255, 200, 200),
                20,
                bottomLeftOrigin=True,
            )

        return im[slicer[0], slicer[1], slicer[2]]


class DirectUndistortedCamera(UndistortedCamera, DirectCamera):
    pass
