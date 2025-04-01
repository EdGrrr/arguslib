from matplotlib import pyplot as plt
import numpy as np
from .instruments import Camera, Position


class CameraArray:
    def __init__(self, cameras: list[Camera], layout_shape: tuple[int, int]):
        self.cameras = cameras
        self.layout_shape = layout_shape

        if self.layout_shape[0] * self.layout_shape[1] < len(self.cameras):
            raise ValueError(
                "The specified layout doesn't have enough spaces to have a full layout."
            )

        self.positions = self.infer_positions()

    def infer_positions(self) -> list[tuple[int, int]]:
        """
        Infer the positions of the cameras in the layout, based on the lat/lon under the camera Camera.Position properties.
        """
        latitudes = [c.position.lat for c in self.cameras]
        longitudes = [c.position.lon for c in self.cameras]

        # Get the min and max lat/lon
        min_lat, max_lat = min(latitudes), max(latitudes)
        min_lon, max_lon = min(longitudes), max(longitudes)

        # get positions for the layout
        positions = []
        for i in range(self.layout_shape[0]):
            for j in range(self.layout_shape[1]):
                positions.append(
                    Position(
                        (
                            np.mean(longitudes)
                            if self.layout_shape[0] == 1
                            else min_lon
                            + (max_lon - min_lon) * i / (self.layout_shape[0] - 1)
                        ),
                        (
                            np.mean(latitudes)
                            if self.layout_shape[1] == 1
                            else min_lat
                            + (max_lat - min_lat) * j / (self.layout_shape[1] - 1)
                        ),
                        0.1,
                    )
                )

        # get a list of distances between each camera and each position
        distances_to_each_pos_for_each_camera = []
        for camera in self.cameras:
            distances_to_each_pos_for_each_camera.append(
                [
                    (np.array(camera.position.target_xyz(pos)) ** 2).sum()
                    for pos in positions
                ]
            )
        distances_to_each_pos_for_each_camera = np.array(
            distances_to_each_pos_for_each_camera
        )

        closest_pos_index = np.argmin(distances_to_each_pos_for_each_camera, axis=1)
        # get the multiply assigned positions
        multiply_assigned = np.unique(
            closest_pos_index[np.unique(closest_pos_index, return_counts=True)[1] > 1]
        )
        if multiply_assigned.size != 0:
            raise ValueError(
                "There are multiple cameras assigned to the same position. Disambiguation not yet implemented."  # TODO:
            )

        # for each camera, return the position within the layout
        camera_positions = []
        for i in range(len(self.cameras)):
            camera_positions.append(
                (
                    (closest_pos_index[i] // self.layout_shape[1]).item(),
                    (closest_pos_index[i] % self.layout_shape[1]).item(),
                )
            )
        return camera_positions

    def show(self, dt, ax=None, replace_ax=None, label_cameras=True, **kwargs):
        """
        Show the camera array on a map.
        """
        if replace_ax is not None:
            # we need to get the positon of the axis in the figure, and replace it with the new array of subplots...
            fig = replace_ax.figure
            pos = replace_ax.get_subplotspec()
            replace_ax.remove()
            subfig = fig.add_subfigure(
                pos,
            )
            axes = subfig.subplots(
                self.layout_shape[0],
                self.layout_shape[1],
                subplot_kw={"projection": "polar"},
            )
        elif ax is None:
            fig, axes = plt.subplots(
                self.layout_shape[0],
                self.layout_shape[1],
                subplot_kw={"projection": "polar"},
                figsize=(8, 8),
                constrained_layout=True,
            )
        else:
            return self.show(dt, replace_ax=ax, label_cameras=label_cameras)

        for i in range(self.layout_shape[0]):
            for j in range(self.layout_shape[1]):
                ax = axes[self.layout_shape[1] - j - 1, i]
                camera = [
                    c
                    for i_cam, c in enumerate(self.cameras)
                    if self.positions[i_cam] == (i, j)
                ]
                if len(camera) == 0:
                    ax.remove()
                    continue
                camera = camera[0]
                camera.show(
                    dt,
                    pos=f"{self.layout_shape[0]}{self.layout_shape[1]}{3*i+j}",
                    ax=ax,
                    **kwargs,
                )
                ax.set_thetagrids(np.arange(0, 360, 45), labels=[])
                if label_cameras:
                    ax.text(
                        0.15,
                        0.85,
                        f"{camera.attrs['camstr']}".replace("-", "--"),
                        va="bottom",
                        ha="right",
                        transform=ax.transAxes,
                        fontsize="small",
                    )

        return axes

    def annotate_positions(self, positions, ax, *args, **kwargs):
        """
        Annotate the positions of the cameras on the map.
        """
        for i in range(self.layout_shape[0]):
            for j in range(self.layout_shape[1]):
                ax_cam = ax[self.layout_shape[1] - j - 1, i]
                camera = [
                    c
                    for i_cam, c in enumerate(self.cameras)
                    if self.positions[i_cam] == (i, j)
                ]
                if len(camera) == 0:
                    continue
                camera = camera[0]
                camera.annotate_positions(positions, ax_cam, *args, **kwargs)
