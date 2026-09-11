# %%
import datetime as dt
from arguslib.camera.camera_array import CameraArray
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
acft = "a320e5"
icaos = [acft]  # only load the aircraft of interest
# icaos = None  # all aircraft

centre_contrail = True

campaign_name = "COBALT"
tlen = 60 * 60

intersection_time_chunk_size = 10  # seconds

radar = Radar.from_config(campaign_name)
cam = UndistortedCamera.from_config("COBALT", "3-7")
cam2 = UndistortedCamera.from_config("COBALT", "5-2")
multicam = CameraArray([cam, cam2], (2, 1))
ai = AutomaticADSBAircraftInterface(radar, winds="era5", wind_filter=10, icaos=icaos)
ri = RadarInterface(radar, cam)

time1 = dt.datetime.fromisoformat("2025-05-01T07:25:06")
time2 = dt.datetime.fromisoformat("2025-05-01T07:30:23")

fig = plt.figure(figsize=(4, 6))
ax = plt.subplot(321, projection="polar")
ax2 = plt.subplot(322, projection="polar")
ai.camera = multicam
# ai.show(time1, ax=ax)
cam_axes = ai.show(
    time2,
    ax=(ax, ax2),
    replace_ax=None,
)

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
ax2.legend(*ax.get_legend_handles_labels(), loc="lower right")
ax.text(0.02, -0.02, time2.isoformat(), transform=ax.transAxes, va="top")


ax_r1 = plt.subplot(312)
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
)
ax_r2 = plt.subplot(313)
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


# # --- Wind shear and stretching -- with some input from claude
# # nb. it seems like both scans were done in the same direction, so there's not
# # any gain from squeezing the scans to "match" the advected geometry.
# from arguslib.instruments import Position
# from arguslib.misc.geo import xy_offset_to_ll, ft_to_km
# import matplotlib as mpl

# SHEAR_LEVELS = ["300hPa", "350hPa"]

# shear_winds = {}
# for lvl in SHEAR_LEVELS:
#     u, v = ERA5WindData(level=lvl, res="0.25grid", linear_interp="time").get_data(
#         np.array([radar.position.lon]),
#         np.array([radar.position.lat]),
#         time1,
#         simple=False,
#     )
#     alt = hPa_to_km(float(lvl.removesuffix("hPa")))
#     shear_winds[alt] = (float(u[0]), float(v[0]))


# def shear_trails(fleet, acft, dtime, shear_winds, tlen=tlen):
#     """One aircraft's unadvected track, re-advected at each level."""
#     i = fleet.get_ids().index(acft)
#     track = fleet.get_tracks_arr(dtime, tlen=tlen, include_time=True)[i]
#     lon, lat, ages = track[:, 0], track[:, 1], track[:, -1]
#     out = {}
#     for alt, (u, v) in shear_winds.items():
#         dlon, dlat = xy_offset_to_ll(lon, lat, u * ages / 1000, v * ages / 1000)
#         out[alt] = ([Position(a, b, alt) for a, b in zip(dlon, dlat)], ages)
#     return out

# norm = mpl.colors.Normalize(min(shear_winds), max(shear_winds))


# def shear_intersections(
#     fleet,
#     acft,
#     radar,
#     t_scan,
#     ax_r,
#     shear_winds,
#     chunk_size=intersection_time_chunk_size,
#     tlen=tlen,
# ):
#     start, end = radar.get_scan_time_bounds(t_scan)

#     step = chunk_size / 2
#     mids = np.arange(start.timestamp() + step, end.timestamp(), step)
#     edges = np.arange(start.timestamp(), end.timestamp() + step, step)

#     for tm, ti, tf in zip(mids, edges[:-2], edges[2:]):
#         trails = shear_trails(
#             fleet, acft, dt.datetime.fromtimestamp(tm), shear_winds, tlen=tlen
#         )
#         for alt, (positions, ages) in trails.items():
#             got = radar.annotate_intersections(
#                 positions,
#                 ages,
#                 t_scan,
#                 ax_r,
#                 time_bounds=(
#                     dt.datetime.fromtimestamp(ti),
#                     dt.datetime.fromtimestamp(tf),
#                 ),
#                 color="r",
#                 marker="+",
#                 s=30,
#             )


# for t, ax_r in ((time1, ax_r1), (time2, ax_r2)):
#     shear_intersections(ai.fleet, acft, radar, t, ax_r, shear_winds)


# alt_true = ft_to_km(
#     ai.fleet.get_data(time1, "alt_geom", tlen=tlen)[acft]["alt_geom"][-1]
# )


def scan_stretch(radar, t_scan, wind, alt_km=8.0, xlim=None, t_ref=None, min_elev=20.0):
    """Linear x-correction taking radar-relative km to air-relative km."""
    r = radar.data_loader.get_pyart_radar(t_scan)
    t = np.asarray(r.time["data"], float)
    t = t - t[0]
    e = np.asarray(r.elevation["data"], float)
    az = float(r.azimuth["data"][0])
    flip = 90 < az < 270

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

    for tm, ti, tf in zip(mids, edges[:-2], edges[2:]):
        positions, ages = ai.get_trail_positions(
            dt.datetime.fromtimestamp(tm), tlen=tlen
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


def panel_limits(radar, t_scan, acft, wind, width_air_km=15.0, alt_km=8.0):
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
    ax_r.set(xlim=xlim, ylim=YLIM)
    ax_r.set_aspect(1 / abs(k))
#     # ax_r.xaxis.set_major_formatter(FuncFormatter(lambda x, _, f=to_air: f"{f(x):.0f}"))

# --- Sensitivity to ERA5 wind error: re-find the crossings with each trail
# point's wind perturbed by `wind_error_frac` of its own speed, in
# `n_wind_error` directions relative to the wind (a ring of errors).
wind_error_frac = 0.05  # None to skip
n_wind_error = 8
if wind_error_frac:
    z = 1 + wind_error_frac * np.exp(
        1j * np.linspace(0, 2 * np.pi, n_wind_error, endpoint=False)
    )
    ring = [
        dict(wind_scale=s, wind_rotate_deg=r)
        for s, r in zip(np.abs(z), np.rad2deg(np.angle(z)))
    ]

    def crossing_xy(t, xlim, **wind_error):
        """(x, y) km of acft's trail crossing the scan at t, or NaNs if none."""
        crossings = ai.find_intersections(
            t,
            xlim=xlim,
            tlen=tlen,
            chunk_size=intersection_time_chunk_size,
            **wind_error,
        ).get(acft)
        if not crossings:
            return np.nan, np.nan
        positions, ages, time_bounds = crossings[0]
        point = radar.get_intersections(
            positions, ages, t, time_bounds=time_bounds, xlim=xlim
        )[0]
        return point["x"], point["y"]

    # Radar: error bar spanning where the ring's crossings land along the scan.
    for t, ax_r in ((time1, ax_r1), (time2, ax_r2)):
        x0, y0 = crossing_xy(t, ax_r.get_xlim())
        xs = [x0] + [crossing_xy(t, ax_r.get_xlim(), **w)[0] for w in ring]
        err = ax_r.errorbar(
            x0,
            y0,
            xerr=[[x0 - np.nanmin(xs)], [np.nanmax(xs) - x0]],
            fmt="none",
            ecolor="r",
            elinewidth=1,
            capsize=3,
        )
    ax_r1.legend(
        leg.legend_handles + [err],
        [text.get_text() for text in leg.get_texts()]
        # A bare % starts a comment under usetex and truncates the label.
        + [
            f"{wind_error_frac * 100:g}"
            + (r"\%" if plt.rcParams["text.usetex"] else "%")
            + " wind error range"
        ],
        loc="lower right",
    )

    # Camera: thin lines for the two ring members displaced furthest to
    # either side of the trail.
    def cross_track_km(base, other):
        """Mean signed offset (km, +ve to the left) of trail `other` from `base`."""
        lon, lat = np.array([[p.lon, p.lat] for p in base]).T
        lon2, lat2 = np.array([[p.lon, p.lat] for p in other]).T
        kx = 111.111 * np.cos(np.deg2rad(lat))
        tx, ty = np.gradient(lon) * kx, np.gradient(lat) * 111.111
        dx, dy = (lon2 - lon) * kx, (lat2 - lat) * 111.111
        return np.nanmean((tx * dy - ty * dx) / np.hypot(tx, ty))

    ai.camera = multicam

    def trail(**wind_error):
        return ai.get_trail_positions(time2, tlen=tlen, **wind_error)[acft][0]

    base = trail()
    offsets = [cross_track_km(base, trail(**w)) for w in ring]
    for k in (np.nanargmin(offsets), np.nanargmax(offsets)):
        ai.plot_trails(
            time2,
            cam_axes,
            tlen=tlen,
            color_icao=False,
            plot_kwargs={"linewidth": 0.5},
            **ring[k],
        )
    ai.camera = radar

plt.gcf().set_size_inches(4, 6)
plt.savefig("distrail_example.jpg", dpi=400, bbox_inches="tight")
plt.show()
# %%
