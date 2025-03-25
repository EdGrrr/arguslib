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
