from arguslib.instruments.instruments import Position
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

class TimestampedFigure(Figure):
    def __init__(self, *args, **kwargs):
        timestamp = kwargs.pop("timestamp", None)
        if timestamp is None:
            raise ValueError("Timestamped figures need a timestamp property.")
        self.timestamp = timestamp
        Figure.__init__(self, *args, **kwargs)
        


# def plot_range_rings(camera, dt, ranges=[10, 20, 30], alt=10, ax=None, **kwargs):
#     # if ax is None:
#     #     ax = plt.gca()

#     range_out = {}
#     kwargs = {"c": "orange", "lw": 0.7} | kwargs
#     for rd in ranges:
#         # get the positions at the same altitude, but at distance rd
#         global_elevs = np.zeros(37)
#         global_azis = np.arange(0, 361, 10)
#         global_dists = np.ones(37) * rd

#         # convert to lat lon positions
#         positions = [
#             camera.position.ead_to_lla(global_elev, global_azi, global_dist)
#             for global_elev, global_azi, global_dist in zip(
#                 global_elevs, global_azis, global_dists
#             )
#         ]

#         # then move them to the correct altitude
#         positions = [Position(pos.lon, pos.lat, alt) for pos in positions]

#         camera.annotate_positions(
#             positions,
#             None,
#             ax=ax,
#             **kwargs,
#         )
#     return range_out

def plot_range_rings(camera, dt, ranges=[10, 20, 30], alt=10, ax=None, **kwargs):
    # This helper function calculates a destination lat/lon using spherical geometry
    def calculate_destination_point(start_lon, start_lat, bearing_deg, distance_km):
        """
        Calculates a destination point given a starting point, bearing,
        and distance along a great-circle path.
        """
        R = 6371.0  # Average Earth radius in km
        
        lat1_rad = np.deg2rad(start_lat)
        lon1_rad = np.deg2rad(start_lon)
        bearing_rad = np.deg2rad(bearing_deg)
        
        d_div_R = distance_km / R
        
        lat2_rad = np.arcsin(np.sin(lat1_rad) * np.cos(d_div_R) +
                         np.cos(lat1_rad) * np.sin(d_div_R) * np.cos(bearing_rad))
        
        lon2_rad = lon1_rad + np.arctan2(np.sin(bearing_rad) * np.sin(d_div_R) * np.cos(lat1_rad),
                                     np.cos(d_div_R) - np.sin(lat1_rad) * np.sin(lat2_rad))
                                     
        return np.rad2deg(lon2_rad), np.rad2deg(lat2_rad)

    range_out = {}
    kwargs = {"c": "orange", "lw": 0.7} | kwargs
    for rd in ranges:
        # Generate points for a ring at a specific ground distance 'rd'
        azimuths_deg = np.arange(0, 361, 10)
        positions = []
        for azi in azimuths_deg:
            # Calculate the destination lat/lon using the correct spherical math
            target_lon, target_lat = calculate_destination_point(
                camera.position.lon, camera.position.lat, azi, rd
            )
            # Create a Position object at the desired altitude
            positions.append(Position(target_lon, target_lat, alt))

        # Annotate the generated positions on the camera image
        camera.annotate_positions(
            positions,
            dt,
            ax=ax,
            **kwargs,
        )
    return range_out


def plot_beam(
    plotting_instrument, radar, elev_azi, dt=None, ax=None, markers=False, **kwargs
):
    if ax is None:
        ax = plt.gca()

    radar_elev, radar_azimuth = elev_azi
    dists = np.logspace(-2, np.log10(15), 100)
    radar_beam_positions = radar.beam(radar_elev, radar_azimuth, dists)

    kwargs = (
        dict(
            markersize=1,
            color="limegreen",
            lw=0.7,
            label=f"elev={radar_elev:.1f}\\textdegree, az={radar_azimuth:.1f}\\textdegree",
            zorder=10,
        )
        | kwargs
    )

    plotting_instrument.annotate_positions(
        radar_beam_positions[:, 0],
        dt=dt,
        ax=ax,
        **kwargs,
    )

    if markers:
        dists = [0, 2, 5, 10, 15]
        for d in dists:
            pt = radar.beam(radar_elev, radar_azimuth, [d])
            plotting_instrument.annotate_positions(pt[:, 0], dt, ax, "ro", markersize=2)

            plotting_instrument.annotate_positions(
                pt[:, 0],
                dt,
                ax,
                f"---{d:.1f} km",
                fontsize=4,
                color="red",
                plotting_method="text",
            )

    return ax


def get_pixel_transform(camera, ax, lr_flip=True):
    from matplotlib.transforms import Affine2D

    # img_size_px = 3040 * camera.scale_factor
    principal_point = camera.intrinsic.principal_point

    translation_px = (
        camera.image_size_px[0] / 2 - principal_point[0],
        camera.image_size_px[1] / 2 - principal_point[1],
    )

    transPixel = (
        Affine2D()
        .translate(*translation_px)
        .scale(1 / camera.image_size_px[0], 1 / camera.image_size_px[1])
        .rotate_deg_around(0.5, 0.5, -1 * camera.rotation[-1])
    )

    try:
        if (
            not lr_flip
        ):  # can't figure out why this needs doing to gett the unflipped version.
            # seems to be that the default is to flip it for polar plots??
            transPixel = transPixel + Affine2D().scale(-1, 1).translate(1, 0)
        elif (ax.get_theta_direction() == np.pi / 2) and (
            ax.get_theta_direction() == -1
        ):
            # bearing axes, so should have been flipped
            print * ("Warning: bearing axes require flipped projection to be accurate")
    except AttributeError:
        # non-polar axes
        pass

    transPixel = transPixel + ax.transAxes
    return transPixel


def make_camera_axes(
    camera, theta_behaviour="bearing", fig=None, pos=111, replace_ax=None, dt=None
):
    if fig is None:
        fig = plt.figure(FigureClass=TimestampedFigure, timestamp=dt)

    if camera.camera_type == "perspective":
        projection = None
    elif camera.camera_type == "allsky":
        projection = "polar"

    if replace_ax is not None:
        fig = replace_ax.figure
    
        # 1. Get the subplot's grid specification (its "slot" in the layout)
        spec = replace_ax.get_subplotspec()
        
        # 2. Remove the old axes from the figure
        replace_ax.remove()
        
        # 3. Add a new subplot in the SAME grid slot
        ax = fig.add_subplot(spec, projection=projection)
    else:
        ax = fig.add_subplot(pos, projection=projection)

    if projection == "polar":
        if theta_behaviour == "pixels":
            ax.set_theta_offset(np.deg2rad(camera.rotation[-1]))
        elif theta_behaviour == "bearing":
            ax.set_theta_offset(np.pi / 2)
            ax.set_theta_direction(-1)
        elif theta_behaviour == "unflipped_ordinal_aligned":
            ax.set_theta_offset(np.pi / 2)
            # ax.set_theta_direction(-1)
        else:
            raise ValueError(
                "theta_behaviour must be one of 'pixels', 'bearing', or 'unflipped_ordinal_aligned'"
            )
    else:
        ax.set_box_aspect(camera.image_size_px[1] / camera.image_size_px[0])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    return ax
