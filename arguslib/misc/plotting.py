import numpy as np
import matplotlib.pyplot as plt


def plot_range_rings(camera, ranges=[10, 20, 30], alt=10, ax=None, **kwargs):
    if ax is None:
        ax = plt.gca()

    range_out = {}
    kwargs = {"c": "orange", "lw": 0.7} | kwargs
    for rd in ranges:
        range_out[rd] = {}
        rl = []
        for az in range(0, 361, 10):
            elev, dist = np.rad2deg(np.arctan2(alt, rd)), np.sqrt(alt**2 + rd**2)
            rl.append(camera.iead_to_pix(*camera.gead_to_iead(elev, az, dist)))
        rl = np.array(rl)
        ax.plot(rl[:, 0], rl[:, 1], **kwargs)
        range_out[rd]["px"] = rl[:, 0]
        range_out[rd]["py"] = rl[:, 1]
    return range_out


def plot_beam(
    plotting_instrument, radar, elev_azi, dt=None, ax=None, markers=False, **kwargs
):
    if ax is None:
        ax = plt.gca()

    radar_elev, radar_azimuth = elev_azi
    dists = np.logspace(-2, np.log10(15), 100)
    radar_beam_positions = radar.beam(radar_elev, radar_azimuth, dists)

    plotting_instrument.annotate_positions(
        radar_beam_positions[:, 0],
        dt=dt,
        ax=ax,
        markersize=1,
        color="limegreen",
        lw=0.7,
        label=f"elev={radar_elev:.1f}\\textdegree, az={radar_azimuth:.1f}\\textdegree",
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

    if ax.name != "polar":
        return ax.transData

    img_size_px = 3040 * camera.scale_factor
    principal_point = camera.intrinsic.principal_point

    translation_px = 3040 / 2 - principal_point[0], 3040 / 2 - principal_point[1]

    transPixel = (
        Affine2D()
        .translate(*translation_px)
        .scale(1 / img_size_px, 1 / img_size_px)
        .rotate_deg_around(0.5, 0.5, -1 * camera.rotation)
    )

    if (
        not lr_flip
    ):  # can't figure out why this needs doing to gett the unflipped version.
        # seems to be that the default is to flip it for polar plots??
        transPixel = transPixel + Affine2D().scale(-1, 1).translate(1, 0)
    elif (ax.get_theta_direction() == np.pi / 2) and (ax.get_theta_direction() == -1):
        # bearing axes, so should have been flipped
        raise print * (
            "Warning: bearing axes require flipped projection to be accurate"
        )

    transPixel = transPixel + ax.transAxes
    return transPixel


def make_camera_axes(
    camera, theta_behaviour="bearing", fig=None, pos=111, replace_ax=None
):
    if fig is None:
        fig = plt.figure()

    if replace_ax is not None:
        fig = replace_ax.figure
        pos = replace_ax.get_position()
        replace_ax.remove()

        ax = fig.add_axes(
            pos,
            projection="polar",
        )
    else:
        ax = fig.add_subplot(pos, projection="polar")

    if theta_behaviour == "pixels":
        ax.set_theta_offset(np.deg2rad(camera.rotation))
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
    return ax
