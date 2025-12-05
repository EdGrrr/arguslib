# %%
import datetime as dt
from arguslib import (
    UndistortedCamera,
    Radar,
    AutomaticADSBAircraftInterface,
    RadarInterface,
)

cam = UndistortedCamera.from_config("COBALT", "3-7")
rad = Radar.from_config("COBALT")

ai = AutomaticADSBAircraftInterface(RadarInterface(rad, cam))
ax = ai.show(
    dt.datetime(2025, 5, 1, 7, 41, 18),
    color_icao=True,
    trail_kwargs=dict(
        label_acft=True
    ),
    tlen=30 * 60,
)
ax[-1].legend()

# ax = ai.show(dt.datetime(2025, 5, 1,7,30,23), color_icao=True, trail_kwargs=dict(plot_kwargs=dict(plotting_method='intersect_plot'), label_acft=True), tlen=30*60)
# ax[-1].legend()
# %%
