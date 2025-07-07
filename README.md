# Argus Library

## Camera Configuration
A config can be used to pre-define some camera calibration details relating them to campaigns, then loaded using the `.from_config()` method on the relevant object class.

Configs are read from `$HOME/.config/arguslib/c`, `$HOME/.arguslib/`, and `/etc/arguslib/`.

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
  2-11:
    focal_lengths: [1712.4324324324323, 856.2162162162161]
    principal_point: [2304, 1296]
    position: [-0.1791071, 51.499189, 0.1]
    rotation: [45, 250, 190]

```
A null value for the calibration file will cause the code to fall back on the default config (which is distributed with the package files).

The hard-coded pattern used to find camera files is:
```
/disk1/Data/ARGUS/{campaign}/{camstr}/videos/{year}-{mon:0>2}-{day:0>2}/argus-{camstr}_{year}{mon:0>2}{day:0>2}_{hour:0>2}{min:0>2}{second:0>2}_A.mp4
```


### `camera_arrays.yml`
Camera arrays have a name, a campaign, a list of cameras, and a layout, which describes the size of the 2-D grid they live on.
```yaml
COBALTArray:
  campaign: "COBALT"
  cameras:
    - 3-7
    - 5-1
    - 5-2
    - 5-3
    - 5-4
  layout_shape: [3, 3]
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




## Camera Calibration
Cameras need intrinsic and extrinsic calibration to be able to geolocate features.

- Intrinsic calibration links the positions in an image to locations relative to the imager. The required parameters are the coefficients of two polynomials to link the distances, and the principal point of the image:
  - `poly_incident_angle_to_radius`
  - `poly_radius_to_z`
  - `principal_point`
- Extrinsic calibration links the camera coordinate system to the global coordinates. For this, three angles are needed. We determine the direction of the image optical axis using the instrument elevation and azimuth. Then, the roll determines the orientation of the image about this axis.
  - `elevation` is measured between the horizontal plane and the optical axis.
  - `azimuth` is measured between the north axis and the projection of the optical axis into the focal plane.
  - The `roll` is measured between the top of the image and the plane defined by the zenith and focal axes. `rotation` is an alias for `roll`.