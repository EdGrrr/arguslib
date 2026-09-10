# %%
import datetime as dt
from arguslib.misc.geo import hPa_to_km
from csat2.ECMWF import ERA5Data
from arguslib.camera.undistorted_camera import UndistortedCamera
from arguslib.radar.radar import Radar
from arguslib.radar.radar_interface import RadarInterface
from arguslib.aircraft import AutomaticADSBAircraftInterface
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

# %%
icao_include = ["a320e5"]  # only show aircraft of interest
# icao_include = None  # Show all aircraft

centre_contrail = True

campaign_name = "COBALT"
tlen = 60 * 60

intersection_time_chunk_size = 5  # seconds

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

wind = tuple(
    float(w)
    for w in ERA5WindData(
        level="350hPa", res="0.25grid", linear_interp="time"
    ).get_data(
        np.array([radar.position.lon]),
        np.array([radar.position.lat]),
        time1,
        simple=False,
    )
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


ax_r1 = plt.subplot(413)
ai.camera = radar
ai.show(
    time1,
    ax=ax_r1,
    var="DBZ",
    colorbar_flag=False,
    title_flag=False,
    tlen=tlen,
    vmax=0,
    vmin=-40,
    trail_kwargs={
        "icao_include": icao_include,
        "intersection_kwargs": {"chunk_size": intersection_time_chunk_size},
    },
)
ax_r1.set(xlim=(-6, 0) if centre_contrail else (-8, 8), ylim=(6.5, 10.5))
ax_r1.set_ylabel("Altitude")
ax_r1.set_xlabel("")
time1_start, time1_end = radar.get_scan_time_bounds(time1)
time1_str = f"{time1_start.isoformat()} to {time1_end.strftime('%H:%M:%S')}"
ax_r1.text(0.02, 0.02, time1_str, transform=ax_r1.transAxes)
leg = ax_r1.legend(
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
ax_r2 = plt.subplot(414)
ai.show(
    time2,
    ax=ax_r2,
    var="DBZ",
    colorbar_flag=False,
    title_flag=False,
    tlen=tlen,
    vmax=0,
    vmin=-40,
    trail_kwargs={
        "icao_include": icao_include,
        "intersection_kwargs": {"chunk_size": intersection_time_chunk_size},
    },
)
ax_r2.set(
    xlim=(-4, 6) if centre_contrail else (-8, 8),
    ylim=(6.5, 10.5),
)
ax_r2.set_ylabel("Altitude")
ax_r2.set_xlabel("")
time2_start, time2_end = radar.get_scan_time_bounds(time2)
time2_str = f"{time2_start.isoformat()} to {time2_end.strftime('%H:%M:%S')}"
ax_r2.text(0.02, 0.02, time2_str, transform=ax_r2.transAxes)


# --- Wind shear and stretching -- with some input from claude
# nb. it seems like both scans were done in the same direction, so there's not
# any gain from squeezing the scans to "match" the advected geometry.
from arguslib.instruments import Position
from arguslib.misc.geo import xy_offset_to_ll, ft_to_km
import matplotlib as mpl

SHEAR_LEVELS = ["300hPa", "350hPa"]

shear_winds = {}
for lvl in SHEAR_LEVELS:
    u, v = ERA5WindData(level=lvl, res="0.25grid", linear_interp="time").get_data(
        np.array([radar.position.lon]),
        np.array([radar.position.lat]),
        time1,
        simple=False,
    )
    alt = hPa_to_km(float(lvl.removesuffix("hPa")))
    shear_winds[alt] = (float(u[0]), float(v[0]))

for alt, (u, v) in sorted(shear_winds.items()):
    print(f"{alt:5.2f} km  u={u:6.1f}  v={v:6.1f}  |w|={np.hypot(u, v):5.1f} m/s")


def shear_trails(fleet, acft, dtime, shear_winds, tlen=tlen):
    """One aircraft's unadvected track, re-advected at each level."""
    i = fleet.get_ids().index(acft)
    track = fleet.get_tracks_arr(dtime, tlen=tlen, include_time=True)[i]
    lon, lat, ages = track[:, 0], track[:, 1], track[:, -1]
    out = {}
    for alt, (u, v) in shear_winds.items():
        dlon, dlat = xy_offset_to_ll(lon, lat, u * ages / 1000, v * ages / 1000)
        out[alt] = ([Position(a, b, alt) for a, b in zip(dlon, dlat)], ages)
    return out


norm = mpl.colors.Normalize(min(shear_winds), max(shear_winds))
acft = icao_include[0]


def shear_intersections(
    fleet,
    acft,
    radar,
    t_scan,
    ax_r,
    shear_winds,
    chunk_size=intersection_time_chunk_size,
    tlen=tlen,
):
    start, end = radar.get_scan_time_bounds(t_scan)
    print(
        f"{t_scan:%H:%M:%S} sweep: {start:%H:%M:%S} to {end:%H:%M:%S} "
        f"({(end - start).total_seconds():.0f} s)"
    )

    step = chunk_size / 2
    mids = np.arange(start.timestamp() + step, end.timestamp(), step)
    edges = np.arange(start.timestamp(), end.timestamp() + step, step)

    for tm, ti, tf in zip(mids, edges[:-2], edges[2:]):
        trails = shear_trails(
            fleet, acft, dt.datetime.fromtimestamp(tm), shear_winds, tlen=tlen
        )
        for alt, (positions, ages) in trails.items():
            got = radar.annotate_intersections(
                positions,
                ages,
                t_scan,
                ax_r,
                time_bounds=(
                    dt.datetime.fromtimestamp(ti),
                    dt.datetime.fromtimestamp(tf),
                ),
                color="r",
                marker="+",
                s=30,
            )

            if got:
                print(
                    f"  {dt.datetime.fromtimestamp(tm):%H:%M:%S} "
                    f"{alt:.1f} km -> {len(got)} crossing(s)"
                )


for t, ax_r in ((time1, ax_r1), (time2, ax_r2)):
    shear_intersections(ai.fleet, acft, radar, t, ax_r, shear_winds)


alt_true = ft_to_km(
    ai.fleet.get_data(time1, "alt_geom", tlen=tlen)[acft]["alt_geom"][-1]
)
print(f"aircraft at {alt_true:.2f} km")


def scan_stretch(radar, t_scan, wind, alt_km=8.0, xlim=None, t_ref=None, min_elev=20.0):
    """Linear x-correction taking radar-relative km to air-relative km."""
    r = radar.data_loader.get_pyart_radar(t_scan)
    print(f"    elev {r.elevation['data'][0]:.1f} -> {r.elevation['data'][-1]:.1f}")
    t = np.asarray(r.time["data"], float)
    t = t - t[0]
    e = np.asarray(r.elevation["data"], float)
    az = float(r.azimuth["data"][0])
    flip = 90 < az < 270
    print(f"    az {az:.1f}  flip={flip}")

    h = alt_km - radar.position.alt
    with np.errstate(divide="ignore"):
        x = h / np.tan(np.deg2rad(e))
    if flip:
        x = -x

    good = np.isfinite(x) & (e >= min_elev)
    if xlim is not None:
        good &= (x > xlim[0]) & (x < xlim[1])
    s, x0 = np.polyfit(t[good], x[good], 1)
    u, v = wind
    u_par = (u * np.sin(np.deg2rad(az)) + v * np.cos(np.deg2rad(az))) / 1000
    if flip:
        u_par = -u_par

    k = 1 - u_par / s
    start, _ = radar.get_scan_time_bounds(t_scan)
    dt_ref = (t_ref - start).total_seconds() if t_ref else 0.0
    c = u_par * (x0 / s + dt_ref)
    return k, (lambda xx: k * xx + c), s, u_par


for t, ax_r, xlim in ((time1, ax_r1, (-8, 2)), (time2, ax_r2, (-4, 6))):
    k, to_air, s, u_par = scan_stretch(radar, t, wind, alt_km=8.0, xlim=xlim)
    print(f"{t:%H:%M:%S}  s={s*1000:6.1f} m/s  u_par={u_par*1000:6.1f} m/s  k={k:5.2f}")
    ax_r.set_aspect(1 / abs(k))
    # ax_r.xaxis.set_major_formatter(FuncFormatter(lambda x, _, f=to_air: f"{f(x):.0f}"))

from arguslib.misc.interpolation import interpolate_to_intersection


def crossing_x_aircraft(
    radar, fleet, acft, t_scan, chunk_size=intersection_time_chunk_size, tlen=tlen
):
    """x (km, radar-relative) where the aircraft's advected trail crosses the sweep."""
    fig, ax_tmp = plt.subplots()
    ax_tmp.set_xlim(-50, 50)  # annotate_intersections filters on ax xlims

    start, end = radar.get_scan_time_bounds(t_scan)
    step = chunk_size / 2
    mids = np.arange(start.timestamp() + step, end.timestamp(), step)
    edges = np.arange(start.timestamp(), end.timestamp() + step, step)
    az = float(radar.data_loader.get_pyart_radar(t_scan).azimuth["data"][0])
    times = radar.data_loader.get_pyart_radar(t_scan).time["data"]
    first_time, last_time = times[0], times[-1]
    print(
        f"...first_time={dt.datetime.fromtimestamp(first_time):%H:%M:%S} "
        f"...last_time={dt.datetime.fromtimestamp(last_time):%H:%M:%S} "
        f"...az={az:.1f}"
        f"...elev={radar.data_loader.get_pyart_radar(t_scan).elevation['data'][0]:.1f}--{radar.data_loader.get_pyart_radar(t_scan).elevation['data'][-1]:.1f}"
    )

    for tm, ti, tf in zip(mids, edges[:-2], edges[2:]):
        positions, ages = ai.get_trail_positions(
            dt.datetime.fromtimestamp(tm), icao_include=[acft], tlen=tlen
        )[acft]
        got = radar.annotate_intersections(
            positions,
            ages,
            t_scan,
            ax_tmp,
            time_bounds=(dt.datetime.fromtimestamp(ti), dt.datetime.fromtimestamp(tf)),
        )
        if got:
            pts = got[0] if isinstance(got[0], (list, tuple, np.ndarray)) else [got[0]]
            eads = np.atleast_2d(radar.position.target_ead(pts))
            elev, azi, dist = eads[0]
            theta = azi % 360 - az % 360
            x = float(dist * np.cos(np.deg2rad(elev)) * np.cos(np.deg2rad(theta)))
            plt.close(fig)
            return -x if 90 < az < 270 else x

    plt.close(fig)
    raise ValueError(f"no crossing found for {t_scan:%H:%M:%S}")


def panel_limits(radar, t_scan, acft, wind, width_air_km=10.0, alt_km=8.0):
    xc = crossing_x_aircraft(radar, ai.fleet, acft, t_scan)
    k, *_ = scan_stretch(radar, t_scan, wind, alt_km=alt_km)
    for _ in range(3):
        half = width_air_km / (2 * abs(k))
        k, to_air, s, u_par = scan_stretch(
            radar, t_scan, wind, alt_km=alt_km, xlim=(xc - half, xc + half)
        )
    half = width_air_km / (2 * abs(k))
    return (xc - half, xc + half), k, to_air, s, u_par


WIDTH_AIR_KM = 10.0
YLIM = (6.5, 10.5)
for t, ax_r in ((time1, ax_r1), (time2, ax_r2)):
    xlim, k, to_air, s, u_par = panel_limits(radar, t, acft, wind, WIDTH_AIR_KM)
    print(
        f"{t:%H:%M:%S}  centre={np.mean(xlim):5.1f} km  "
        f"s={s*1000:6.1f}  u_par={u_par*1000:6.1f}  k={k:5.2f}  "
        f"xlim=({xlim[0]:.1f}, {xlim[1]:.1f})"
    )
    ax_r.set(xlim=xlim, ylim=YLIM)
    ax_r.set_aspect(1 / abs(k))
    # ax_r.xaxis.set_major_formatter(FuncFormatter(lambda x, _, f=to_air: f"{f(x):.0f}"))

plt.gcf().set_size_inches(4, 7)
plt.savefig("distrail_example.jpg", dpi=200, bbox_inches="tight")
plt.show()
# %%
