# %%
import datetime as dt
from arguslib.camera.undistorted_camera import UndistortedCamera
from arguslib.radar.radar import Radar
from arguslib.radar.radar_interface import RadarInterface
from arguslib.aircraft import AutomaticADSBAircraftInterface
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

icao_include = ["a320e5"]  # only show aircraft of interest
# icao_include = None  # Show all aircraft

campaign_name = "COBALT"
tlen = 60 * 60

radar = Radar.from_config(campaign_name)
cam = UndistortedCamera.from_config("COBALT", "3-7")
ai = AutomaticADSBAircraftInterface(radar, winds="era5", wind_filter=10)
ri = RadarInterface(radar, cam)

time1 = dt.datetime.fromisoformat("2025-05-01T07:25:06")
time2 = dt.datetime.fromisoformat("2025-05-01T07:30:23")

ax = plt.subplot(211, projection="polar")
ai.camera = cam
# ai.show(time1, ax=ax)
ai.show(time2, ax=ax, trail_kwargs={"icao_include": icao_include})

from csat2.ECMWF import ERA5WindData

wind = ERA5WindData(level="200hPa", res="0.25grid").get_data(
    radar.position.lon, radar.position.lat, time1
)

ri.annotate_sampled_scan(time2, ax, wind=(0, 0), color="grey", label="Scan pattern")
ri.annotate_sampled_scan(
    time1,
    ax,
    wind=wind,
    t_ref=time2,
    color="darkgreen",
    label="Scan 1 (advection-corrected)",
)
ri.annotate_sampled_scan(
    time2, ax, wind=wind, color="limegreen", label="Scan 2 (no advection correction)"
)
ax.legend(loc="lower right")
ax.text(0.02, 0.02, time2.isoformat(), transform=ax.transAxes)


ax = plt.subplot(413)
ai.camera = radar
ai.show(
    time1,
    ax=ax,
    var="DBZ",
    colorbar_flag=False,
    title_flag=False,
    tlen=tlen,
    vmax=0,
    vmin=-40,
    trail_kwargs={"icao_include": icao_include},
)
ax.set(xlim=(-9, 9), ylim=(6, 12))
ax.set_ylabel("Altitude")
ax.set_xlabel("")
time1_start, time1_end = radar.get_scan_time_bounds(time1)
time1_str = f"{time1_start.isoformat()} to {time1_end.strftime('%H:%M:%S')}"
ax.text(0.02, 0.02, time1_str, transform=ax.transAxes)
leg = ax.legend(
    [
        Line2D(
            [0],
            [0],
            marker="x",
            color="r",
            label="Aircraft position (advection-corrected)",
            markersize=5,
            linestyle="None",
            linewidth=0,
        )
    ],
    ["Aircraft position (advection-corrected)"],
    loc="upper right",
)
ax = plt.subplot(414)
ai.show(
    time2,
    ax=ax,
    var="DBZ",
    colorbar_flag=False,
    title_flag=False,
    tlen=tlen,
    vmax=0,
    vmin=-40,
    trail_kwargs={"icao_include": icao_include},
)
ax.set(xlim=(-9, 9), ylim=(6, 12))
ax.set_ylabel("Altitude")
ax.set_xlabel("")
time2_start, time2_end = radar.get_scan_time_bounds(time2)
time2_str = f"{time2_start.isoformat()} to {time2_end.strftime('%H:%M:%S')}"
ax.text(0.02, 0.02, time2_str, transform=ax.transAxes)

plt.gcf().set_size_inches(4, 7)
plt.savefig("distrail_example.jpg", dpi=200, bbox_inches="tight")
plt.show()


# %%
