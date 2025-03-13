from calibration import Projection, unit
from geo import haversine, bearing
import numpy as np


class Position:
    # Position as
    # - lon, lat, alt (assuming Earth is locally flat)
    def __init__(self, lon, lat, alt):
        self.lon = lon
        self.lat = lat
        self.alt = alt

    # Elevation, Azimuth, distance
    #  - Elevation from local horizontal
    #  - Azimuth relative to local north (compass bearing)
    #  - Distance in km

    def target_ead(self, target_position):
        # Approximate - can use pyproj for exact (or with large distances)
        distance = haversine(
            self.lon, self.lat, target_position.lon, target_position.lat
        )
        alt_diff = target_position.alt - self.alt
        target_distance = np.sqrt(distance**2 + alt_diff**2)
        target_elevation = np.rad2deg(np.arctan2(alt_diff, distance))
        target_azimuth = bearing(
            self.lon, self.lat, target_position.lon, target_position.lat
        )
        return np.array([target_elevation, target_azimuth, target_distance])

    # XYZ
    #  - X East-West distace in km
    #  - Y North-south distance in km
    #  - Z altitude in km

    def target_xyz(self, target_position):
        # Assume Earth is locally flat, but also spherical
        distance = haversine(
            self.lon, self.lat, target_position.lon, target_position.lat
        )
        target_z = target_position.alt - self.alt
        target_bearing = bearing(
            self.lon, self.lat, target_position.lon, target_position.lat
        )
        target_x = distance * np.sin(np.deg2rad(target_bearing))
        target_y = distance * np.cos(np.deg2rad(target_bearing))
        return np.array([target_x, target_y, target_z])

    def ead_to_xyz(self, target_elevation, target_azimuth, target_distance):
        target_z = target_distance * np.sin(np.deg2rad(target_elevation))
        target_hor_range = target_distance * np.cos(np.deg2rad(target_elevation))
        target_x = target_hor_range * np.sin(np.deg2rad(target_azimuth))
        target_y = target_hor_range * np.cos(np.deg2rad(target_azimuth))
        return np.array([target_x, target_y, target_z])

    def ead_to_lla(self, target_elevation, target_azimuth, target_distance):
        return self.xyz_to_lla(
            *self.ead_to_xyz(target_elevation, target_azimuth, target_distance)
        )

    def xyz_to_ead(self, target_x, target_y, target_z):
        alt_diff = target_position.alt - self.alt
        target_distance = np.sqrt(target_x**2 + target_y**2 + target_z**2)
        target_elevation = np.rad2deg(
            np.arctan2(target_z, np.sqrt(target_x**2 + target_y**2))
        )
        target_azimuth = np.rad2deg(np.arctan2(target_x, target_y))
        return np.array([target_elevation, target_azimuth, target_distance])

    def xyz_to_lla(self, target_x, target_y, target_z):
        lon = self.lon + target_x / (111.111 * np.cos(self.lat))
        lat = self.lat + target_y / (111.111)
        return Position(lon, lat, target_z + self.alt)

    def __repr__(self):
        return f"{self.lon:.3f} {self.lat:.3f} {self.alt:.3f}"


class Instrument:
    def __init__(self, position: Position, rotation: list):
        self.position = position
        self.rotation = rotation

    # lla and xyz are functions of the instrument position, rather than the instrument itself
    # ead on the other hand is an instrument property, as it depends on the instrument rotation.

    def iead_to_gead(self, elevation, azimuth, dist):
        """Instrument El, az, dist to global values, based on an
        instrument facing directly up/north.  For the allsky cameras,
        this is probably more straightforwrd. A little trickier for
        perspective cameras.

        """
        return elevation, azimuth + self.rotation, dist

    def gead_to_iead(self, elevation, azimuth, dist):
        return elevation, azimuth - self.rotation, dist

    def target_iead(self, target_position: Position):
        ead = self.position.target_ead(target_position)
        return self.gead_to_iead(*ead)

    def iead_to_lla(self, elevation, azimuth, dist):
        return self.position.ead_to_lla(*self.iead_to_gead(elevation, azimuth, dist))

    def iead_to_xyzworld(self, elevation, azimuth, dist):
        return self.position.ead_to_xyz(*self.iead_to_gead(elevation, azimuth, dist))

    def xyzworld_to_iead(target_x, target_y, target_z):
        gead = self.position.xyz_to_ead(target_x, target_y, target_z)
        return self.gead_to_iead(*gead)


class Camera(Instrument):
    def __init__(
        self,
        intrinsic_calibration: Projection,
        *args,
        scale_factor=1,
        camera_type="allsky",
        **kwargs,
    ):
        self.intrinsic = intrinsic_calibration
        self.scale_factor = scale_factor
        self.camera_type = camera_type
        super().__init__(*args, **kwargs)

    @classmethod
    def from_filename(cls, filename, *args, **kwargs):
        return cls(Projection.fromfile(filename), *args, **kwargs)

    # view here is in the camera space - we need to get the xyz in camera projection, not the world projection
    def target_pix(self, target_position: Position):
        return self.iead_to_pix(*self.target_iead(target_position))

    def pix_to_iead(self, pix_x, pix_y, distance):
        xyz = self.intrinsic.image_to_view(
            [pix_x * self.scale_factor, pix_y * self.scale_factor]
        )
        # xyz for projection, convert to ead for projection
        return self.position.xyz_to_ead(distance * unit(xyz))

    def iead_to_pix(self, elevation, azimuth, dist=10):
        return (
            self.intrinsic.view_to_image(
                self.position.ead_to_xyz(elevation, azimuth, dist)
            )
            / self.scale_factor
        )

    def radar_beam(self, target_elevation, target_azimuth, radar):
        output = []
        for beam_r in range(0, 361, 30):
            output.append(
                self.el_azi_to_pixel(
                    target_elevation + radar.beamwidth * np.cos(np.deg2rad(beam_r)),
                    target_azimuth + radar.beamwidth * np.sin(np.deg2rad(beam_r)),
                    self.ppx,
                    self.ppy,
                )
            )
        return np.array(output)


class Radar(Instrument):
    def __init__(self, beamwidth, *args, **kwargs):
        self.beamwidth = beamwidth
        super().__init__(*args, **kwargs)

    def beam(self, radar_elevation, radar_azimuth, radar_distances):
        if self.beamwidth:
            return np.array(
                [
                    [
                        self.iead_to_lla(radar_elevation, radar_azimuth, radar_dist),
                        self.iead_to_lla(
                            radar_elevation + self.beamwidth / 2,
                            radar_azimuth,
                            radar_dist,
                        ),
                        self.iead_to_lla(
                            radar_elevation,
                            radar_azimuth + self.beamwidth / 2,
                            radar_dist,
                        ),
                        self.iead_to_lla(
                            radar_elevation - self.beamwidth / 2,
                            radar_azimuth,
                            radar_dist,
                        ),
                        self.iead_to_lla(
                            radar_elevation,
                            radar_azimuth - self.beamwidth / 2,
                            radar_dist,
                        ),
                    ]
                    for radar_dist in radar_distances
                ]
            )
        else:
            return np.array(
                [
                    [self.iead_to_lla(radar_elevation, radar_azimuth, radar_dist)]
                    for radar_dist in radar_distances
                ]
            )
