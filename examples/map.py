# %%
import datetime as dt

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from arguslib import (
    MapInstrument,
    CameraArray,
    AutomaticADSBAircraftInterface,
    UndistortedCamera,
    Position,
)
from arguslib.misc.plotting import plot_range_rings

import matplotlib.pyplot as plt

# 1. Create the map instrument to serve as the background
map_instrument = MapInstrument.from_config("chilbolton_and_sirta")
ai = AutomaticADSBAircraftInterface(map_instrument)
ax = ai.show(
    dt.datetime(2025, 5, 1, 7, 25, 6),
    tlen=20 * 60,
    trail_kwargs={"plot_kwargs": {"linewidth": 0.5}},
)

# 2. Load an instrument to use as the center for the rings
multicam = CameraArray.from_config("COBALTArray", camera_class=UndistortedCamera)

paris_cam = UndistortedCamera.from_config("COBALT", "3-7")
paris_cam.position = Position(2.3514, 48.8575, 0.1)


# 3. Plot the range rings on the map
#    - plotting_instrument is the map
#    - center_instrument is the camera
plot_range_rings(
    plotting_instrument=map_instrument,
    center_instrument=multicam,
    dt=None,
    ax=ax,
    ranges=[10],
    label="10 km range (COBALT)",
)
plot_range_rings(
    plotting_instrument=map_instrument,
    center_instrument=paris_cam,
    dt=None,
    ax=ax,
    ranges=[10],
    label="10 km range (Hypothetical Paris Camera)",
    c="teal",
)


ax.scatter(
    [-0.17899],
    [51.49911],
    s=4,
    color="blue",
    label="Huxley Building, Imperial College London",
    transform=ccrs.PlateCarree(),
    zorder=9,
)
ax.scatter(
    [2.35591],
    [48.84718],
    s=4,
    color="navy",
    label="Sorbonne Sciences, Paris",
    transform=ccrs.PlateCarree(),
    zorder=9,
)

plt.legend(fontsize="small")

# %%
ai.camera = multicam
ai.show(dt.datetime(2025, 5, 1, 7, 25, 6), tlen=20 * 60)
# %%
ax = map_instrument.show(dt.datetime(2025, 5, 1, 7, 25, 6))
plot_range_rings(
    plotting_instrument=map_instrument,
    center_instrument=multicam,
    dt=None,
    ax=ax,
    ranges=[10],
    label="10 km range",
)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)

multicam.show(dt.datetime(2025, 5, 1, 7, 25, 6))

# %%
