import math
import numpy as np

R_E = 6372.8  # see http://rosettacode.org/wiki/Haversine_formula

def haversine(lon1, lat1, lon2, lat2):
    """
    Computes the Haversine distance between two points in km.
    This version is vectorized to handle numpy arrays.
    """
    # Ensure inputs are numpy arrays for vector operations
    lon1, lat1, lon2, lat2 = map(np.asarray, [lon1, lat1, lon2, lat2])

    dlat = np.deg2rad(lat2 - lat1)
    dlon = np.deg2rad(lon2 - lon1)

    lat1_rad = np.deg2rad(lat1)
    lat2_rad = np.deg2rad(lat2)

    arc = np.sin(dlat / 2) ** 2 + np.cos(lat1_rad) * np.cos(lat2_rad) * (np.sin(dlon / 2) ** 2)
    c = 2 * np.arcsin(np.sqrt(arc))
    return R_E * c

def bearing(lon1, lat1, lon2, lat2):
    """Calculates the bearing between two points
    NOTE: switched lon/lat order to match TASIC code

    All arguments in degrees

    Code from: https://www.igismap.com/formula-to-find-bearing-or-heading-angle-between-two-points-latitude-longitude/"""
    dlon = np.deg2rad(lon2 - lon1)
    X = np.cos(np.deg2rad(lat2)) * np.sin(dlon)
    Y = np.cos(np.deg2rad(lat1)) * np.sin(np.deg2rad(lat2)) - np.sin(
        np.deg2rad(lat1)
    ) * np.cos(np.deg2rad(lat2)) * np.cos(dlon)
    return np.rad2deg(np.arctan2(X, Y))

def xy_offset_to_ll(lon1, lat1, xoff, yoff):
    return (
        lon1 + (xoff / (111.111 * np.cos(np.deg2rad(lat1)))),
        lat1 + (yoff / 111.111),
    )
    
def destination_point(lon, lat, bearing_deg, distance_km):
    """
    Calculates a destination point given a starting point, bearing,
    and distance along a great-circle path on a sphere.
    """
    R = 6371.0  # Average Earth radius in km
    
    lat1_rad = np.deg2rad(lat)
    lon1_rad = np.deg2rad(lon)
    bearing_rad = np.deg2rad(bearing_deg)
    
    d_div_R = distance_km / R
    
    # Clamp the argument to arcsin to avoid NaN from floating point errors
    sin_lat1 = np.sin(lat1_rad)
    cos_lat1 = np.cos(lat1_rad)
    sin_d_R = np.sin(d_div_R)
    cos_d_R = np.cos(d_div_R)
    
    arg_for_arcsin = sin_lat1 * cos_d_R + cos_lat1 * sin_d_R * np.cos(bearing_rad)
    lat2_rad = np.arcsin(np.clip(arg_for_arcsin, -1.0, 1.0))
    
    lon2_rad = lon1_rad + np.arctan2(np.sin(bearing_rad) * sin_d_R * cos_lat1,
                                 cos_d_R - sin_lat1 * np.sin(lat2_rad))
                                 
    return np.rad2deg(lon2_rad), np.rad2deg(lat2_rad)


def calculate_bearing(lat1, lon1, lat2, lon2):
    # Convert from degrees to radians
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)

    dlon = lon2 - lon1

    x = np.sin(dlon) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)

    initial_bearing = np.arctan2(x, y)
    
    # Convert to degrees and normalize
    bearing = (np.degrees(initial_bearing) + 360) % 360
    return bearing
    

def hPa_to_ft(press):
    return 145366.45 * (1 - (press / 1013.25) ** 0.190284)  # 100x the original

def ft_to_km(ft):
    return ft * 0.0003048

def hPa_to_km(press):
    return ft_to_km(hPa_to_ft(press))

def ft_to_hPa(ft):
    meters = ft * 0.3048
    return 1013.25 * (1 - (meters / 44330.0)) ** 5.255

def km_to_hPa(km):
    meters = km * 1000
    return 1013.25 * (1 - (meters / 44330.0)) ** 5.255


def xr_add_cyclic_points(da):
    """
    Add cyclic points at start and end of `lon` dimension of data array.
    
    From StackExchange.
    
    Inputs
    da: xr.DataArray including dimensions (lat,lon)
    """
    import xarray as xr
    # Borrows heavily from cartopy.util.add_cyclic point, but adds at start and end.

    lon_idx = da.dims.index('lon')
    
    start_slice = [slice(None)] * da.ndim
    end_slice = [slice(None)] * da.ndim
    start_slice[lon_idx] = slice(0, 1)
    end_slice[lon_idx] = slice(-1, None)
    
    wrap_data = np.concatenate([da.values[tuple(end_slice)], da.values, da.values[tuple(start_slice)]], axis=lon_idx)
    wrap_lon = np.concatenate([da.lon.values[-1:] - 360, da.lon.values, da.lon.values[0:1] + 360])

    # Generate output DataArray with new data but same structure as input
    outp_da = xr.DataArray(data=wrap_data, 
                           coords=dict(lat=da.lat, lon=wrap_lon), 
                           dims=da.dims, 
                           attrs=da.attrs)
    
    return outp_da