.. _configuration:

Configuration
=============

`arguslib` uses a set of configuration files to define instruments, data paths, and other parameters. This allows you to adapt the library to different campaigns and setups without changing the code.

Configuration File Locations
----------------------------

Configuration files are read from the following locations, in order. Settings from files found later in the list will override settings from earlier ones.

1.  `/etc/arguslib/`
2.  `$HOME/.arguslib/`
3.  `$HOME/.config/arguslib/`


Configuration File Reference
----------------------------

Here are the primary configuration files used by `arguslib`.

### `cameras.yml`

This file defines the properties of individual camera instruments. Cameras are grouped by a `campaign` name and identified by a unique `camstr`.

There are two types of cameras you can configure: **all-sky (fisheye)** and **perspective**.

**Common Properties:**

*   `position`: A list of `[longitude, latitude, altitude_km]` for the camera's location.

**All-sky Camera Properties:**

*   `calibration_file`: Path to a YAML file containing the camera's intrinsic polynomial calibration. If set to `null`, a default calibration file distributed with the package is used.
*   `rotation`: A single number representing the roll of the camera in degrees, which orients the image relative to North.

**Perspective Camera Properties:**

*   `focal_lengths`: A list `[fx, fy]` defining the focal length in pixels.
*   `principal_point`: A list `[cx, cy]` defining the image center in pixels.
*   `distortion_coeffs` (optional): A list of distortion coefficients `[k1, k2, p1, p2, k3]`.
*   `rotation`: A list of three angles `[elevation, azimuth, roll]` in degrees that define the camera's pointing direction and orientation.

**Example `cameras.yml`:**

```yaml
COBALT:
  # Example of an all-sky (fisheye) camera
  3-7:
    calibration_file:  ~/.config/arguslib/cam_3-7_calibration.yml
    position: [-1.439252, 51.146668, 0.1]
    rotation: 10 # Single roll angle

  # Example of a perspective camera
  2-11:
    focal_lengths: [1712.43, 856.21]
    principal_point: [2304, 1296]
    position: [-0.1791071, 51.499189, 0.1]
    rotation: [45, 250, 190] # [elevation, azimuth, roll]
```

### `camera_arrays.yml`

This file groups multiple cameras (defined in `cameras.yml`) into a single named array for simultaneous visualization.

*   `campaign`: The campaign the cameras belong to.
*   `cameras`: A list of `camstr` identifiers to include in the array.
*   `layout_shape`: A tuple `[rows, cols]` defining the grid for displaying the cameras.

**Example `camera_arrays.yml`:**

```yaml
COBALT_Array:
  campaign: "COBALT"
  cameras:
    - "3-7"
    - "5-1"
    - "5-2"
    - "5-3"
    - "5-4"
  layout_shape: [3, 3]
```

### `radars.yml`

This file defines the properties of the radar for a given campaign.

*   `beamwidth`: The radar's beamwidth in degrees.
*   `position`: A list of `[longitude, latitude, altitude_km]`.
*   `rotation`: A single number for the radar's rotation/boresight alignment.

**Example `radars.yml`:**

```yaml
COBALT:
  beamwidth: 0.6
  position: [-1.439339, 51.144752, 0.1]
  rotation: 0
```

### `adsb_path.txt`

This is a simple text file used by the `AircraftInterface` to find ADS-B flight data. It should contain a single line with the absolute path to the directory where the data files are stored.

**Example `adsb_path.txt`:**

```
/disk1/Data/ADS-B/COBALT/
```

Data File Locations
-------------------

The library currently expects data files to be organized in specific directory structures. The paths are constructed using hard-coded patterns.

*   **Camera Video Files:**
    `/disk1/Data/ARGUS/{campaign}/{camstr}/videos/{year}-{mon:0>2}-{day:0>2}/argus-{camstr}_{year}{mon:0>2}{day:0>2}_{hour:0>2}{min:0>2}{second:0>2}_A.mp4`
*   **Radar RHI Files:**
    `/disk1/Data/{campaign}/radar/L1/{year}{mon:0>2}{day:0>2}/ncas-mobile-ka-band-radar-1_cao_{year}{mon:0>2}{day:0>2}-{hour:0>2}**_rhi_l1_v1.0.0.nc`
*   **Radar VPT Files:** (Limited support)
    `/disk1/Data/{campaign}/radar/L1/{year}{mon:0>2}{day:0>2}/ncas-mobile-ka-band-radar-1_cao_{year}{mon:0>2}{day:0>2}-{hour:0>2}**_vpt_l1_v1.0.0.nc`
*   **ADS-B Flight Data:**
    The files `YYYYMMDD_ADS-B.nc` and `YYYYMMDD_ADS-B.txt` are expected to be in the directory specified in `adsb_path.txt`.


Camera Calibration Concepts
---------------------------

Cameras need intrinsic and extrinsic calibration to be able to geolocate features.

- Intrinsic calibration links the positions in an image to locations relative to the imager. The required parameters are the coefficients of two polynomials to link the distances, and the principal point of the image:
  - `poly_incident_angle_to_radius`
  - `poly_radius_to_z`
  - `principal_point`
- Extrinsic calibration links the camera coordinate system to the global coordinates. For this, three angles are needed. We determine the direction of the image optical axis using the instrument elevation and azimuth. Then, the roll determines the orientation of the image about this axis.
  - `elevation` is measured between the horizontal plane and the optical axis.
  - `azimuth` is measured between the north axis and the projection of the optical axis into the focal plane.
  - The `roll` is measured between the top of the image and the plane defined by the zenith and focal axes. `rotation` is an alias for `roll`.