# Argus Library

## Camera Configuration
A config can be used to pre-define some camera calibration details relating them to campaigns.
Configs are read from `$HOME/.config/arguslib/c`, `$HOME/.arguslib/`, and `/etc/arguslib/`

Some example config files are below.

In the current implementation, the same campaign and instrument identifying string are required in the pattern matching for loading data.


### `cameras.yml`
Cameras belong to a campaign and are identified using a `camstr`.
```yaml
COBALT:
  3-7:
    calibration_file:  ~/.arguslib/cam1_calibration.yml
    position: [-1.439252, 51.146668, 0.1]
    rotation: 10
  3-8:
    calibration_file:  null
    position: [-1.438511, 51.149064, 0.1]
    rotation: 160
  5-1:
    calibration_file:  null
    position: [-1.363346, 50.952369, 0.1]
    rotation: 70
  5-2:
    calibration_file:  null
    position: [-1.141453, 51.147445, 0.1]
    rotation: -80
  5-3:
    calibration_file:  null
    position: [-1.443658, 51.320375, 0.1]
    rotation: 45
  5-4:
    calibration_file:  null
    position: [-1.78568, 51.09967, 0.1]
    rotation: -115

```
A null value for the calibration file will cause the code to fall back on the default config (which is distributed with the package files).

The hard-coded pattern used to find camera files is:
```
/disk1/Data/ARGUS/{campaign}/{camstr}/videos/{year}-{mon:0>2}-{day:0>2}/argus-{camstr}_{year}{mon:0>2}{day:0>2}_{hour:0>2}{min:0>2}{second:0>2}_A.mp4
```


### `cameras.yml`
Cameras belong to a campaign and are identified using a `camstr`.
```yaml
COBALT:
  3-7:
    calibration_file:  ~/.arguslib/cam1_calibration.yml
    position: [-1.439252, 51.146668, 0.1]
    rotation: 15
  3-8:
    calibration_file:  null
    position: [-1.438511, 51.149064, 0.1]
    rotation: 160
  5-1:
    calibration_file:  null
    position: [-1.363346, 50.952369, 0.1]
    rotation: 70
  5-2:
    calibration_file:  null
    position: [-1.141453, 51.147445, 0.1]
    rotation: -80
  5-3:
    calibration_file:  null
    position: [-1.443658, 51.320375, 0.1]
    rotation: 45
  5-4:
    calibration_file:  null
    position: [-1.78568, 51.09967, 0.1]
    rotation: -115

```
A null value for the calibration file will cause the code to fall back on the default config (which is distributed with the package files).

The hard-coded pattern used to find camera MP4 files containing the image data is:
```sh
/disk1/Data/ARGUS/{campaign}/{camstr}/videos/{year}-{mon:0>2}-{day:0>2}/argus-{camstr}_{year}{mon:0>2}{day:0>2}_{hour:0>2}{min:0>2}{second:0>2}_A.mp4
```

### `radars.yml`
Currently, one radar belongs to each campaign.
```yaml
COBALT:
  beamwidth: 0.6
  position: [-1.439339, 51.144752, 0.1]
  rotation: 0
``` 

The hard-coded patterns used to find processed radar `vpt` and `rhi` netcdfs are:
```sh
/disk1/Data/{campaign}/radar/L1/{year}{mon:0>2}{day:0>2}/ncas-mobile-ka-band-radar-1_cao_{year}{mon:0>2}{day:0>2}-{hour:0>2}**_vpt_l1_v1.0.0.nc
/disk1/Data/{campaign}/radar/L1/{year}{mon:0>2}{day:0>2}/ncas-mobile-ka-band-radar-1_cao_{year}{mon:0>2}{day:0>2}-{hour:0>2}**_rhi_l1_v1.0.0.nc
```
`vpt` data has limited support.
