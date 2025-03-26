import numpy as np
import matplotlib.pyplot as plt

from arguslib.arguslib.instruments.instruments import Camera


def plot_range_rings(camera: Camera, ranges=[10, 20, 30], alt=10, ax=None):
    if ax is None:
        ax = plt.gca()

    range_out = {}
    for rd in ranges:
        range_out[rd] = {}
        rl = []
        for az in range(0, 361, 10):
            elev, dist = np.rad2deg(np.arctan2(alt, rd)), np.sqrt(alt**2 + rd**2)
            rl.append(camera.iead_to_pix(*camera.gead_to_iead(elev, az, dist)))
        rl = np.array(rl)
        ax.plot(rl[:, 0], rl[:, 1], c="orange", lw=0.7)
        range_out[rd]["px"] = rl[:, 0]
        range_out[rd]["py"] = rl[:, 1]
    return range_out


def plot_beam(camera: Camera, radar, elev_azi, ax=None, **kwargs):
    if ax is None:
        ax = plt.gca()

    radar_elev, radar_azimuth = elev_azi
    dists = np.logspace(-1, 1.0, 100)
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

    return ax
