# Argus Library

## Camera Configuration
A config can be used to pre-define some camera calibration details relating them to campaigns.
Configs are read from `$HOME/.config/arguslib/cameras.yml`, `$HOME/.arguslib/cameras.yml`, and `/etc/arguslib/cameras.yml`

An example config is below.
```yaml
COBALT: # campaign
  3-7: # camstr
    calibration_file:  ~/.arguslib/cam1_calibration.yml
    position: [-1.439252, 51.146668, 0.1]
    rotation: 10
  3-8:
    calibration_file:  null
    position: [-1.439252, 51.146668, 0.1]
    rotation: 10
```
A null value for the calibration file will cause the code to fall back on the default config (which is in the package files).

These can be used for creating `arguslib.instruments.Camera` instances.
In the current implementation, the same campaign and camstr are required in the pattern matching for loading data in `arguslib.video.CameraData` objects.
The hard-coded pattern is:
```
/disk1/Data/ARGUS/{campaign}/{camstr}/videos/{year}-{mon:0>2}-{day:0>2}/argus-{camstr}_{year}{mon:0>2}{day:0>2}_{hour:0>2}{min:0>2}{second:0>2}_A.mp4
```
