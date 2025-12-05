from arguslib.instruments.instruments import Position
from matplotlib.figure import Figure
from ..misc.geo import destination_point

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import transforms
from datetime import timezone  # added


class TimestampedFigure(Figure):
    def __init__(self, *args, **kwargs):
        timestamp = kwargs.pop("timestamp", None)
        if timestamp is None:
            raise ValueError("Timestamped figures need a timestamp property.")
        self.timestamp = timestamp
        Figure.__init__(self, *args, **kwargs)


def plot_range_rings(
    plotting_instrument,
    center_instrument,
    dt,
    ranges=[10, 20, 30],
    alt=10,
    ax=None,
    **kwargs,
):
    """
    Plots range rings on a plottable instrument's axes.

    Each ring is plotted in a separate call to ensure they are not connected.

    Args:
        plotting_instrument (PlottableInstrument): The instrument to plot on (e.g., MapInstrument, Camera).
        center_instrument (Instrument): The instrument providing the center position for the rings.
        dt (datetime): The datetime for the plot.
        ranges (list): List of ranges in km.
        alt (float): Altitude for the ring positions.
        ax (Axes): The axes to plot on.
        **kwargs: Keyword arguments for plotting.
    """

    if hasattr(center_instrument, "cameras"):
        label = kwargs.pop("label", None)

        for i, cam in enumerate(center_instrument.cameras):
            plot_range_rings(
                plotting_instrument,
                cam,
                dt,
                ranges=ranges,
                alt=alt,
                ax=ax,
                label=label if i == 0 else None,
                **kwargs,
            )
        return ax

    plot_kwargs = {"c": "orange", "lw": 0.7} | kwargs
    azimuths_deg = np.arange(0, 361, 10)

    # Loop to generate and plot each ring separately
    for rd in ranges:
        # Vectorized calculation for all points in a single ring
        # target_lons, target_lats = destination_point(
        #     center_instrument.position.lon,
        #     center_instrument.position.lat,
        #     azimuths_deg,
        #     rd,
        # )
        elev = np.rad2deg(np.arctan2(alt, rd))  # in degrees
        dist = np.sqrt(rd**2 + alt**2)

        positions = [
            center_instrument.position.ead_to_lla(elev, azimuth, dist)
            for azimuth in azimuths_deg
        ]

        # A single annotation call for each ring
        plotting_instrument.annotate_positions(
            positions,
            dt,
            ax=ax,
            **plot_kwargs,
        )

    return ax


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
    import numpy as np
    from arguslib.instruments.instruments import rotation_matrix_i_to_g
    from arguslib.camera.calibration import unit

    principal_point = camera.intrinsic.principal_point
    cx, cy = principal_point

    # Determine rotation based on camera type and orientation
    R_i_to_g = rotation_matrix_i_to_g(*camera.rotation)
    # Calculate Global Up in Instrument Frame
    i_up = R_i_to_g.T @ np.array([0, 0, 1])

    # If the camera is looking nearly straight up or down (within ~25 degrees),
    # align to North (bearing). Otherwise, align to Zenith (up is up).
    if np.abs(i_up[2]) > 0.9:
        top_bearing_deg = camera.get_bearing_to_image_top()
        rotation_angle = -top_bearing_deg
    else:
        # Calculate angle of Zenith relative to Image Top (Y-axis)
        # atan2(x, y) gives angle from Y-axis (CW if x is Right)
        zenith_angle = np.degrees(np.arctan2(i_up[0], i_up[1]))
        rotation_angle = 180 - zenith_angle

    w, h = camera.image_size_px

    if camera.camera_type == "perspective":
        axes_aspect = h / w
        scale_factor = 1 / w
    else:
        axes_aspect = 1.0
        scale_factor = 1 / min(w, h)

    transPixel = (
        Affine2D()
        .translate(-cx, -cy)
        .scale(scale_factor, scale_factor)
        .rotate_deg(rotation_angle)
        .scale(1, 1 / axes_aspect)
        .translate(0.5, 0.5)
    )

    try:
        if not lr_flip:
            transPixel = transPixel + Affine2D().scale(-1, 1).translate(1, 0)
    except AttributeError:
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


def get_timestamp_from_ax(ax):
    # ax is an axes on a timestamped figure, or an "axes iterable" which can contain either axes or more axes itereables.
    try:
        # ax is a matplotlib ax
        fig = ax.get_figure()
        if fig is None:
            return None
        return fig.timestamp
    except AttributeError as e:
        if "SubFigure" in e.args[0]:
            return ax.get_figure().get_figure().timestamp
        try:
            i = 0
            timestamp = None
            while timestamp is None and i < len(ax):
                timestamp = get_timestamp_from_ax(ax[i])
                i += 1
            return timestamp
        except TypeError:
            return None
        
def get_fig_from_ax_or_axs(ax):
    # ax is an axes on a timestamped figure, or an "axes iterable" which can contain either axes or more axes itereables.
    try:
        # ax is a matplotlib ax
        fig = ax.get_figure()
        return fig
    except AttributeError as e:
        if "SubFigure" in e.args[0]:
            return ax.get_figure().get_figure()
        try:
            i = 0
            fig = None
            while fig is None and i < len(ax):
                fig = get_fig_from_ax_or_axs(ax[i])
                i += 1
            return fig
        except TypeError:
            return None


# --- Solar position helpers (NOAA/SPA-style approximation) ---
def _wrap_deg(x):
    return (x % 360.0 + 360.0) % 360.0


def _jd_utc(dt):
    # Expect UTC; if tz-aware convert to UTC, else assume already UTC
    if dt.tzinfo is not None:
        dt = dt.astimezone(timezone.utc)
    year = dt.year
    month = dt.month
    day = dt.day
    hour = dt.hour
    minute = dt.minute
    second = dt.second + dt.microsecond / 1e6

    a = (14 - month) // 12
    y = year + 4800 - a
    m = month + 12 * a - 3
    jdn = day + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045
    day_frac = (hour - 12) / 24 + minute / 1440 + second / 86400
    return jdn + day_frac


def _solar_alt_az(dt, lat_deg, lon_deg):
    # Based on NOAA/SPA approximations, good to ~0.01 deg
    jd = _jd_utc(dt)
    T = (jd - 2451545.0) / 36525.0

    L0 = _wrap_deg(280.46646 + 36000.76983 * T + 0.0003032 * T * T)  # mean long
    M = _wrap_deg(357.52911 + 35999.05029 * T - 0.0001537 * T * T)  # mean anomaly
    e = 0.016708634 - 0.000042037 * T - 0.0000001267 * T * T  # eccentricity

    M_rad = np.deg2rad(M)
    C = (
        (1.914602 - 0.004817 * T - 0.000014 * T * T) * np.sin(M_rad)
        + (0.019993 - 0.000101 * T) * np.sin(2 * M_rad)
        + 0.000289 * np.sin(3 * M_rad)
    )
    true_long = L0 + C
    Omega = 125.04 - 1934.136 * T
    lambda_app = true_long - 0.00569 - 0.00478 * np.sin(np.deg2rad(Omega))

    # Obliquity
    eps0 = (
        23.0
        + 26.0 / 60.0
        + 21.448 / 3600.0
        - (46.8150 / 3600.0) * T
        - (0.00059 / 3600.0) * T * T
        + (0.001813 / 3600.0) * T * T * T
    )
    eps = eps0 + 0.00256 * np.cos(np.deg2rad(Omega))

    # RA/Dec
    lam_rad = np.deg2rad(lambda_app)
    eps_rad = np.deg2rad(eps)
    sin_lam = np.sin(lam_rad)
    cos_lam = np.cos(lam_rad)
    sin_eps = np.sin(eps_rad)
    cos_eps = np.cos(eps_rad)

    alpha = np.rad2deg(np.arctan2(cos_eps * sin_lam, cos_lam))  # degrees
    alpha = _wrap_deg(alpha)
    delta = np.rad2deg(np.arcsin(sin_eps * sin_lam))  # degrees

    # Sidereal time (GMST) -> LST
    GMST = _wrap_deg(
        280.46061837
        + 360.98564736629 * (jd - 2451545.0)
        + 0.000387933 * T * T
        - (T**3) / 38710000.0
    )
    LST = _wrap_deg(GMST + lon_deg)

    # Hour angle
    H = _wrap_deg(LST - alpha)
    if H > 180:
        H -= 360.0

    # Convert to alt/az
    lat_rad = np.deg2rad(lat_deg)
    dec_rad = np.deg2rad(delta)
    H_rad = np.deg2rad(H)

    sin_alt = np.sin(lat_rad) * np.sin(dec_rad) + np.cos(lat_rad) * np.cos(
        dec_rad
    ) * np.cos(H_rad)
    alt_rad = np.arcsin(sin_alt)

    # Azimuth from North, clockwise
    y = np.sin(H_rad)
    x = np.cos(H_rad) * np.sin(lat_rad) - np.tan(dec_rad) * np.cos(lat_rad)
    az_rad = np.arctan2(y, x)
    az_deg = _wrap_deg(
        np.rad2deg(az_rad) + 180.0
    )  # shift from South-based to North-based

    alt_deg = np.rad2deg(alt_rad)

    # Optional: simple atmospheric refraction correction (sea level)
    if alt_deg > -1.0:
        R = 1.02 / np.tan(np.deg2rad(alt_deg + 10.3 / (alt_deg + 5.11))) / 60.0
        alt_deg += R

    return alt_deg, az_deg


# --- end helpers ---


def show_the_sun(ax, camera, **kwargs):
    dt = ax.get_figure().timestamp

    lat = camera.position.lat
    lon = camera.position.lon

    sun_elevation, sun_azimuth = _solar_alt_az(dt, lat, lon)

    solar_position_approx = camera.position.ead_to_lla(sun_elevation, sun_azimuth, 1)

    camera.annotate_positions(
        [solar_position_approx],
        dt,
        ax=ax,
        **kwargs,
    )
    return ax
