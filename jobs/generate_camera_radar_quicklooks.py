# %%
import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import tqdm

from arguslib.aircraft import AircraftInterface
from arguslib.radar import RadarInterface
from arguslib.instruments.radar import Radar
from arguslib.instruments.camera_array import CameraArray

adsb_datadir = Path("/disk1/Data/ADS-B/COBALT/")
outdir = Path(__file__).parent / "output" / "camera_radar_quicklooks"
outdir.mkdir(exist_ok=True, parents=True)

multicam = CameraArray.from_config("COBALTArray")
radar = Radar.from_config("COBALT")
cri = RadarInterface(radar, multicam)
cai = AircraftInterface(cri)

start_time = datetime.datetime(
    2025,
    3,
    28,
    0,
    0,
)
dt = start_time


# tqdm object, but dont estimate the remaining because it is not known
pbar = tqdm.tqdm(
    total=None,
    desc="Generating quicklooks",
    unit="quicklook",
)

# loop over all times in the radar data
while True:
    next_radar_time = radar.data_loader.get_next_time(dt)
    if next_radar_time is None:
        break
    else:
        dt = next_radar_time
        # print("Next radar time:", dt)
    try:
        adsb_file = adsb_datadir / (dt.strftime("%Y%m%d") + "_ADS-B")
        if cai.fleet.loaded_file != str(adsb_file):
            cai.fleet.load_output(str(adsb_file))
    except FileNotFoundError:
        # print("  Skipping due to missing ADS-B data.")
        continue
    except RuntimeError as e:
        if "NetCDF: HDF error" in str(e):
            # Some broken ADS-B files, we need to skip to the next day...
            dt = dt + datetime.timedelta(days=1)
            dt.replace(hour=0, minute=0, second=0)
            continue

    try:
        axes_cams, ax_radar = cai.show(
            dt,
            tlen=60 * 60,
            color_icao=True,
            trail_kwargs={
                "plot_kwargs": {"cam_kwargs": {"max_range_km": 30, "alpha": 0.5}}
            },
        )
    except FileNotFoundError as e:
        if "No video file found" in str(e) or "No camera data found" in str(e):
            # expected e.g. during nighttime
            # print("  Skipping due to missing camera data.")
            continue
        else:
            raise e
    except TypeError as e:
        # This seems to happen sometimes when the radar data isn't available.
        # log and continue
        if "see help(pcolormesh)" in str(e):
            # print("  Skipping due to missing radar data.")
            print(f"  Skipping {dt} due to missing radar data.", e)
            continue
        else:
            raise e

    fig = ax_radar.figure

    day_outdir = outdir / dt.strftime("%Y%m%d")
    day_outdir.mkdir(exist_ok=True, parents=True)
    fig.savefig(
        day_outdir / f"quicklook_{dt.strftime('%Y%m%d_%H%M%S')}.png",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.1,
    )
    plt.close(fig)

    pbar.update(1)

print("Done. Final time:", dt)

# %%
