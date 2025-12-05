# %%
from arguslib import DirectUndistortedCamera, VideoInterface, AircraftInterface
import datetime
from pathlib import Path

cam = DirectUndistortedCamera.from_config("COBALT", "3-7")
dt_start = datetime.datetime(2025, 5, 11, 7, 30)
outdir = Path(__file__).parent / "output" / "videos_from_interface"
outdir.mkdir(exist_ok=True, parents=True)

# %%
# --- Using VideoInterface for camera-only video ---
video_iface = VideoInterface(cam)
output_video_file = outdir / f"cam_video_{dt_start.strftime('%Y-%m-%d_%H%M')}.mp4"

num_frames = 90
video_iface.generate_video(
    output_path=str(output_video_file),
    start_dt=dt_start,
    end_dt=dt_start + datetime.timedelta(minutes=num_frames - 1),
    step_timedelta=datetime.timedelta(minutes=1),
    fps=4,  # From your original example
    show_kwargs={"brightness_adjust": 1.0},  # Optional: kwargs for DirectCamera.show
    time_overlay=True,  # Optional: to add timestamp on frames
)

# %%
# --- A video with annotated trails ---
aci = AircraftInterface(cam)  # cam is a DirectUndistortedCamera
aci.load_flight_data(dt_start)
video_iface_with_trails = VideoInterface(aci)  # Uses the same 'cam' instance

output_video_file_trails = (
    outdir / f"trails_video_{dt_start.strftime('%Y-%m-%d_%H%M')}.mp4"
)


num_frames = 90
video_iface_with_trails.generate_video(
    output_path=str(output_video_file_trails),
    start_dt=dt_start,
    end_dt=dt_start + datetime.timedelta(minutes=num_frames - 1),
    step_timedelta=datetime.timedelta(minutes=1),
    fps=4,
    show_kwargs={
        "brightness_adjust": 1.0,
        "tlen": 15 * 60,
    },  # Optional: kwargs for DirectCamera.show
    time_overlay=True,  # Optional: to add timestamp on frames
)
# takes some time to render...

# %%
