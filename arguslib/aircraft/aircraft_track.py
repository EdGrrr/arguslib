# %%
from pathlib import Path
from PIL import Image
import requests
import numpy as np
from io import BytesIO
import matplotlib.pyplot as plt
from arguslib.misc import geo
from arguslib.instruments import Camera, Radar, Position


cam1_calibration_file = (
    Path(__file__).parent.parent / "instruments/cam1_calibration.yml"
)
cam1_position = Position(-1.439252, 51.146668, 0.1)
cam1 = Camera.from_filename(
    cam1_calibration_file.absolute(),
    cam1_position,
    rotation=10,
    scale_factor=3040 / 500,
)

# camra_position = Position(-1.43812, 51.144980, 0.13)
# camra_position_1km = Position(-1.43812, 51.144980, 1.1)
# xband_position = Position(-1.43552, 51.146055, 0.10)

# Get aircraft data
response = requests.get(
    f"https://opendata.adsb.fi/api/v2/lat/{cam1_position.lat}/lon/{cam1_position.lon}/dist/100"
)
aircraft = response.json()

# Get latest image
session = requests.Session()
session.auth = ("argus", "panoptes")
img_response = session.get("http://argus.edgryspeerdt.com/static/argus-3-7.jpg")
img = Image.open(BytesIO(img_response.content))

# img = Image.open('argus-3-7.jpg')


plt.imshow(img)

range_out = {}
# Plot range rings
ranges = [10, 20, 30]
alt = 10
for rd in ranges:
    range_out[rd] = {}
    rl = []
    for az in range(0, 361, 10):
        elev, dist = np.rad2deg(np.arctan2(alt, rd)), np.sqrt(alt**2 + rd**2)
        rl.append(cam1.iead_to_pix(*cam1.gead_to_iead(elev, az, dist)))
    rl = np.array(rl)
    plt.plot(rl[:, 0], rl[:, 1], c="orange", lw=0.7)
    range_out[rd]["px"] = rl[:, 0]
    range_out[rd]["py"] = rl[:, 1]

acft_out = {}

# Plot aircraft positions
for acft in aircraft["aircraft"]:
    try:
        if acft["alt_geom"] < 26000:
            continue
        dist = geo.haversine(
            cam1.position.lon, cam1.position.lat, acft["lon"], acft["lat"]
        )
        acft_out[acft["hex"]] = {
            "acft_dist": dist,
            "flight": acft["flight"],
            "t": acft["t"],
        }

        if dist < 90:
            pl = cam1.target_pix(
                Position(acft["lon"], acft["lat"], acft["alt_geom"] / (3.33 * 1000))
            )
            plt.scatter(pl[0], pl[1], c="r", s=2)
            acft_out[acft["hex"]]["acft_px"] = pl[0]
            acft_out[acft["hex"]]["acft_py"] = pl[1]

        # Plot trails

        # gs is in knots, not 100% about wind, but probably also knots
        wind_u, wind_v = (
            0.51444 * acft["ws"] * np.sin(np.deg2rad(acft["wd"])),
            0.51444 * acft["ws"] * np.cos(np.deg2rad(acft["wd"])),
        )
        acft_u, acft_v = (
            0.51444 * acft["gs"] * np.sin(np.deg2rad(acft["track"])),
            0.51444 * acft["gs"] * np.cos(np.deg2rad(acft["track"])),
        )

        # To get the track location, we wind back the aircraft position,
        # then integrate the wind forwards. A simple solution in then
        times = np.arange(-3600, 0.1, 10)  # Seconds
        track_offset_km = (
            np.array([acft_u - wind_u, acft_v - wind_v])[:, None] * times / 1000
        )
        track_pos = np.array(
            geo.xy_offset_to_ll(acft["lon"], acft["lat"], *track_offset_km)
        )
        pl_track = np.array(
            [
                cam1.target_pix(
                    Position(tp[0], tp[1], acft["alt_geom"] / (3.33 * 1000))
                )
                for tp in track_pos.T
            ]
        )

        dist = geo.haversine(
            cam1.position.lon, cam1.position.lat, track_pos[0], track_pos[1]
        )
        plt.plot(pl_track.T[0][dist < 90], pl_track.T[1][dist < 90], c="r", lw=0.5)

        acft_out[acft["hex"]]["track_px"] = pl_track.T[0][dist < 90]
        acft_out[acft["hex"]]["track_py"] = pl_track.T[1][dist < 90]

    except KeyError:
        pass


# acft_out2 = statuslib.aircraft_trails(cam1)


# Get radar direction
radar_lon, radar_lat = -1.439339, 51.144752
radar_beamwidth = 0.6
kepler = Radar(radar_beamwidth, Position(radar_lon, radar_lat, 0.1), rotation=0)
radar_elev, radar_azimuth = 40, -180

# Plot radar direction
radar_beam_locs = kepler.beam(radar_elev, radar_azimuth, [0, 8, 35])

radar_out = {}

radar_beam_px = np.array([cam1.target_pix(pt) for pt in radar_beam_locs[1]])
plt.plot(radar_beam_px[:, 0], radar_beam_px[:, 1], c="green")
radar_out["radar_near_px"] = radar_beam_px[:, 0]
radar_out["radar_near_py"] = radar_beam_px[:, 1]

radar_beam_px = np.array([cam1.target_pix(pt) for pt in radar_beam_locs[2]])
plt.plot(radar_beam_px[:, 0], radar_beam_px[:, 1], c="green")
radar_out["radar_far_px"] = radar_beam_px[:, 0]
radar_out["radar_far_py"] = radar_beam_px[:, 1]

radar_beam_px = np.array([cam1.target_pix(pt) for pt in radar_beam_locs[:, 0]])
plt.plot(radar_beam_px[:, 0], radar_beam_px[:, 1], c="green")
radar_out["radar_ray_px"] = radar_beam_px[:, 0]
radar_out["radar_ray_py"] = radar_beam_px[:, 1]

plt.show()

# %%
