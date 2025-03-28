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


def plot_beam(camera, radar, elev_azi, ax=None, markers=True, **kwargs):
    if ax is None:
        ax = plt.gca()

    radar_elev, radar_azimuth = elev_azi
    dists = np.logspace(-2, np.log10(15), 100)
    radar_beam_positions = radar.beam(radar_elev, radar_azimuth, dists)
    radar_beam_pix = np.array(
        [camera.target_pix(pt) for pt in radar_beam_positions.reshape(-1)]
    ).reshape(-1, 5, 2)

    kwargs = {
        "label": f"{radar_elev:.1f}\\textdegree~elevation; {radar_azimuth:.1f}\\textdegree~azimuth",
        **kwargs,
    }
    ax.plot(radar_beam_pix[:, 0, 0], radar_beam_pix[:, 0, 1], **kwargs)

    plot_range_rings(camera, ax=ax)

    if markers:
        dists = [0, 2, 5, 10, 15]
        for d in dists:
            pt = radar.beam(radar_elev, radar_azimuth, [d])[0, 0]
            pix = camera.target_pix(pt)
            ax.plot(pix[0], pix[1], "ro", markersize=2)
            ax.text(pix[0] + 30, pix[1], f"\\qquad{d:.1f} km", fontsize=4, color="red")

    return ax


def get_pixel_transform(camera, ax, lr_flip=True):
    from matplotlib.transforms import Affine2D

    if ax.name != "polar":
        return ax.transData

    img_size_px = 3040 * camera.scale_factor
    principal_point = camera.intrinsic.principal_point

    # does this need to be reversed?
    translation_px = 3040 / 2 - principal_point[0], 3040 / 2 - principal_point[1]

    transPixel = (
        Affine2D()
        .translate(*translation_px)
        .scale(1 / img_size_px, 1 / img_size_px)
        .rotate_deg_around(0.5, 0.5, -1 * camera.rotation)
    )

    if lr_flip:
        transPixel = transPixel + Affine2D().scale(-1, 1).translate(1, 0)
    elif (ax.get_theta_direction() == np.pi / 2) and (ax.get_theta_direction() == -1):
        # bearing axes, so should have been flipped
        raise print * (
            "Warning: bearing axes require flipped projection to be accurate"
        )

    transPixel = transPixel + ax.transAxes
    return transPixel


def make_camera_axes(camera, fig=None, pos=111, theta_behaviour="bearing"):
    if fig is None:
        fig = plt.figure()
    ax = fig.add_subplot(pos, projection="polar")

    if theta_behaviour == "pixels":
        ax.set_theta_offset(-1 * np.deg2rad(camera.rotation))
    elif theta_behaviour == "bearing":
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
    return ax


def plot_rhi_beam(ax, elevation, **kwargs):
    # plot the sweep extremes
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()

    # plot start
    endpoint_start_x = xlims[1], xlims[1] * np.tan(np.deg2rad(elevation))
    endpoint_start_y = (
        ylims[1] * np.tan(np.deg2rad(90 - elevation)),
        ylims[1],
    )
    endpoint = (
        endpoint_start_y
        if endpoint_start_x[1] > endpoint_start_y[0]
        else endpoint_start_x
    )
    ax.plot([0, endpoint[0]], [0, endpoint[1]], c="limegreen", lw=0.7)

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
