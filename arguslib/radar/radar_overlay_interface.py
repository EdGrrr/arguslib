import datetime
import numpy as np
import pyart
from pyart.util import datetime_from_radar

from arguslib.instruments.instruments import PlottableInstrument, Position
# We will duck-type the radar, but it's typically arguslib.instruments.radar.Radar
from arguslib.instruments.camera import Camera # Needed for from_campaign helper
from arguslib.instruments.camera_array import CameraArray # Needed for from_campaign helper
from arguslib.instruments.direct_camera import DirectUndistortedCamera # Needed for from_campaign helper
from arguslib.aircraft.aircraft_interface import AircraftInterface # Needed for from_campaign helper

class RadarOverlayInterface(PlottableInstrument):
    """
    An interface to overlay radar-derived information (like beams or scan locations)
    onto another PlottableInstrument's display (e.g., a Camera).
    """

    def __init__(self, radar: PlottableInstrument, target_instrument: PlottableInstrument):
        """
        Args:
            radar: The radar instrument providing beam/scan information.
                   Must have 'beam', 'data_loader', 'position', 
                   'iead_to_lla', 'beamwidth' attributes/methods.
            target_instrument: The instrument onto which the radar information will be overlaid.
                               This instrument must implement show() and annotate_positions().
        """
        # Duck-typing check for radar capabilities
        required_attrs = ['beam', 'data_loader', 'position', 'iead_to_lla', 'beamwidth', 'attrs']
        for attr in required_attrs:
            if not hasattr(radar, attr):
                raise TypeError(f"Radar-like object must have a '{attr}' attribute/method.")
        if not callable(getattr(radar, 'beam')):
            raise TypeError("Radar-like object 'beam' attribute must be callable.")
        if not callable(getattr(radar, 'iead_to_lla')):
            raise TypeError("Radar-like object 'iead_to_lla' attribute must be callable.")

        self.radar = radar
        self.target_instrument = target_instrument

        # Ensure radar data loader is initialized
        if self.radar.data_loader is None:
            if hasattr(self.radar, 'initialise_data_loader') and callable(getattr(self.radar, 'initialise_data_loader')):
                self.radar.initialise_data_loader()
            else:
                raise AttributeError("Radar-like object has no data_loader and no initialise_data_loader method.")

        # Inherit attributes from the target instrument, potentially adding radar info
        attrs = getattr(self.target_instrument, 'attrs', {}).copy()
        attrs['radar'] = getattr(self.radar, 'attrs', {}) # radar.attrs should exist due to check
        super().__init__(**attrs)

    def show(self, dt: datetime.datetime, ax: any = None, **kwargs: any):
        """
        Shows the target instrument's display. Radar overlays are added via
        separate annotation methods.

        Args:
            dt: Datetime for the display.
            ax: Matplotlib axes to plot on (if target_instrument uses them).
                Expected to be None if target_instrument is a DirectCamera.
            **kwargs: Additional arguments passed to target_instrument.show().

        Returns:
            The axes object returned by target_instrument.show() (or None).
        """
        return self.target_instrument.show(dt, ax=ax, **kwargs)

    def annotate_positions(self, positions: list[Position], dt: datetime.datetime, ax: any = None, *args: any, **kwargs: any):
        """
        Annotates arbitrary positions onto the target instrument's display.
        Delegates to the target instrument's annotate_positions method.
        """
        return self.target_instrument.annotate_positions(positions, dt, ax, *args, **kwargs)

    def annotate_radar_beams(self, dt: datetime.datetime, ax: any = None,
                             ranges_km: np.ndarray = None,
                             beam_type: str = 'start_end', 
                             **kwargs):
        """
        Annotates radar beams onto the target instrument's display.

        Args:
            dt: Datetime for the radar scan/beam.
            ax: Matplotlib axes (or None for DirectCamera).
            ranges_km: Array of distances (km) along the beam for annotation points.
                       If None, uses radar's max range or a default.
            beam_type: Type of beam(s) to annotate:
                       'start_end': Annotate the first and last beams of the scan (default).
                       'active': Annotate the single ray closest to 'dt'.
                       'all_sweeps': Annotate the first ray of each sweep.
            **kwargs: Additional arguments passed to target_instrument.annotate_positions().
        """
        radar_pyart_obj = self.radar.data_loader.get_pyart_radar(dt)

        if ranges_km is None:
            if hasattr(radar_pyart_obj, 'range') and 'data' in radar_pyart_obj.range and len(radar_pyart_obj.range['data']) > 0:
                 max_range_m = radar_pyart_obj.range['data'][-1]
                 ranges_km = np.linspace(0.1, max_range_m / 1000., 50)
            else:
                 ranges_km = np.linspace(0.1, 100, 50) 

        beams_to_annotate = [] 

        if not radar_pyart_obj.elevation['data'].size or not radar_pyart_obj.azimuth['data'].size:
            print(f"Warning: Radar data for {dt} has no elevation or azimuth data. Cannot annotate beams.")
            return ax

        if beam_type == 'start_end':
            beams_to_annotate.append((radar_pyart_obj.elevation['data'][0], radar_pyart_obj.azimuth['data'][0]))
            if len(radar_pyart_obj.elevation['data']) > 1:
                 beams_to_annotate.append((radar_pyart_obj.elevation['data'][-1], radar_pyart_obj.azimuth['data'][-1]))
        elif beam_type == 'active':
            if hasattr(radar_pyart_obj, 'time') and 'data' in radar_pyart_obj.time and len(radar_pyart_obj.time['data']) > 0:
                ray_times_seconds = radar_pyart_obj.time['data']
                scan_start_datetime = datetime_from_radar(radar_pyart_obj)
                dt_seconds_since_scan_start = (dt - scan_start_datetime).total_seconds()
                closest_ray_index = np.argmin(np.abs(ray_times_seconds - dt_seconds_since_scan_start))
                if 0 <= closest_ray_index < len(radar_pyart_obj.elevation['data']):
                     beams_to_annotate.append((radar_pyart_obj.elevation['data'][closest_ray_index],
                                               radar_pyart_obj.azimuth['data'][closest_ray_index]))
                else:
                     print(f"Warning: Could not find active ray index {closest_ray_index} for {dt}. No active beam annotated.")
            else:
                 print(f"Warning: Radar object has no time data. Cannot determine active beam for {dt}. No active beam annotated.")
        elif beam_type == 'all_sweeps':
             if hasattr(radar_pyart_obj, 'sweep_start_ray_index') and 'data' in radar_pyart_obj.sweep_start_ray_index:
                 for sweep_index in range(radar_pyart_obj.nsweeps):
                     start_ray_index = radar_pyart_obj.sweep_start_ray_index['data'][sweep_index]
                     if start_ray_index < len(radar_pyart_obj.elevation['data']):
                          beams_to_annotate.append((radar_pyart_obj.elevation['data'][start_ray_index],
                                                    radar_pyart_obj.azimuth['data'][start_ray_index]))
             else:
                  print(f"Warning: Radar object has no sweep information. Cannot annotate all sweeps for {dt}. No beams annotated.")
        else:
            print(f"Warning: Unknown beam_type '{beam_type}'. No beams annotated.")
            return ax

        last_ax = ax
        default_colors = ['darkgreen', 'limegreen', 'blue', 'purple', 'orange']
        for i, (el, az) in enumerate(beams_to_annotate):
            beam_lla_cross_sections = self.radar.beam(el, az, ranges_km)
            centerline_lla_points = [cross_section[0] for cross_section in beam_lla_cross_sections]
            beam_kwargs = kwargs.copy()
            if 'color' not in beam_kwargs and 'c' not in beam_kwargs:
                 beam_kwargs['color'] = default_colors[i % len(default_colors)]
            if 'label' not in beam_kwargs:
                 if beam_type == 'start_end': beam_kwargs['label'] = 'Start Beam' if i == 0 else 'End Beam'
                 elif beam_type == 'active': beam_kwargs['label'] = 'Active Beam'
                 elif beam_type == 'all_sweeps': beam_kwargs['label'] = f'Sweep {i+1} Start'
            last_ax = self.target_instrument.annotate_positions(centerline_lla_points, dt, ax, **beam_kwargs)
        return last_ax

    def annotate_scan_box(self, dt: datetime.datetime, ax: any = None, range_km: float = 10.0, **kwargs):
        """ Annotates a box representing the cross-scan extent. """
        radar_pyart_obj = self.radar.data_loader.get_pyart_radar(dt)
        if not radar_pyart_obj.elevation['data'].size or not radar_pyart_obj.azimuth['data'].size:
            print(f"Warning: No elevation or azimuth data in radar scan for {dt}. Cannot annotate scan box.")
            return ax

        start_el, start_az = radar_pyart_obj.elevation['data'][0], radar_pyart_obj.azimuth['data'][0]
        end_el, end_az = radar_pyart_obj.elevation['data'][-1], radar_pyart_obj.azimuth['data'][-1]
        
        # Use the radar's beamwidth attribute
        if not hasattr(self.radar, 'beamwidth') or self.radar.beamwidth is None:
            print("Warning: Radar has no beamwidth attribute. Cannot calculate scan box accurately.")
            return ax

        orthogonal_direction = (start_az + 90) % 360
        center_start_pos = self.radar.iead_to_lla(start_el, start_az, range_km)
        cross_scan_dist_km = range_km * np.tan(np.deg2rad(self.radar.beamwidth / 2.))
        
        corner1 = center_start_pos.ead_to_lla(0, orthogonal_direction, cross_scan_dist_km)
        corner2 = center_start_pos.ead_to_lla(0, orthogonal_direction, -cross_scan_dist_km)
        center_end_pos = self.radar.iead_to_lla(end_el, end_az, range_km)
        corner3 = center_end_pos.ead_to_lla(0, orthogonal_direction, -cross_scan_dist_km)
        corner4 = center_end_pos.ead_to_lla(0, orthogonal_direction, cross_scan_dist_km)
        
        position_corners = [corner1, corner2, corner3, corner4, corner1]
        box_kwargs = {'color': 'yellow', 'linewidth': 0.2, 'label': f'{range_km:.0f} km Scan Box'}
        box_kwargs.update(kwargs)
        return self.target_instrument.annotate_positions(position_corners, dt, ax=ax, **box_kwargs)

    @staticmethod
    def _create_plottable_instrument_from_config(campaign: str, config_dict: dict):
        target_type = config_dict.get('type')
        camstr = config_dict.get('camstr') # Common for Camera and DirectCamera

        if target_type == 'Camera':
            if camstr is None: raise ValueError("Camera config requires 'camstr'.")
            return Camera.from_config(campaign, camstr)
        elif target_type == 'DirectCamera':
             if camstr is None: raise ValueError("DirectCamera config requires 'camstr'.")
             return DirectUndistortedCamera.from_config(campaign, camstr)
        elif target_type == 'CameraArray':
             array_name = config_dict.get('array_name')
             if array_name is None: raise ValueError("CameraArray config requires 'array_name'.")
             return CameraArray.from_config(array_name)
        elif target_type == 'AircraftInterface':
             camera_config = config_dict.get('camera_config')
             if camera_config is None: raise ValueError("AircraftInterface config requires 'camera_config'.")
             camera_instrument = RadarOverlayInterface._create_plottable_instrument_from_config(campaign, camera_config)
             # AircraftInterface can take a PlottableInstrument, which RadarOverlayInterface is.
             # However, it's more typical for AircraftInterface to wrap a base camera.
             # If the intention is to overlay radar on an aircraft view, the AircraftInterface should wrap the camera,
             # and then *that* AircraftInterface becomes the target_instrument for RadarOverlayInterface.
             # For this factory, we assume camera_config defines the base camera for AircraftInterface.
             return AircraftInterface(camera_instrument)
        else:
            raise ValueError(f"Unsupported target_instrument type: {target_type}")

    @classmethod
    def from_campaign(cls, campaign: str, target_instrument_config: dict):
        """
        Factory method to create RadarOverlayInterface.
        The radar is always the default radar for the campaign.
        The target_instrument is created based on target_instrument_config.
        """
        from arguslib.instruments.radar import Radar # Import here to avoid circular dependency at module level
        radar_instrument = Radar.from_config(campaign)
        target_instrument = cls._create_plottable_instrument_from_config(campaign, target_instrument_config)
        return cls(radar_instrument, target_instrument)

    def to_image_array(self, time=True):
        """
        If the target_instrument supports to_image_array (e.g., is a DirectCamera),
        this method calls its to_image_array() method.
        """
        if hasattr(self.target_instrument, 'to_image_array') and callable(getattr(self.target_instrument, 'to_image_array')):
            return self.target_instrument.to_image_array(time=time)
        else:
            raise NotImplementedError(
                "to_image_array is only available if the target_instrument supports it (e.g., is a DirectCamera or wraps one)."
            )

    @property
    def image(self):
        if hasattr(self.target_instrument, 'image'):
            return self.target_instrument.image
        else:
            raise NotImplementedError(
                "image property is only available if the target_instrument supports it (e.g., is a DirectCamera or wraps one)."
            )