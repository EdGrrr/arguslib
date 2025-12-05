# %%
import arguslib
import arguslib.misc
import arguslib.misc.thermo
import datetime
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

variables = [
    "lon",
    "lat",
    "alt_baro",
    "alt_geom",
    "geom_rate",
    "tas",
    "gs",
    "ws",
    "track",
    "true_heading",
    "wd",
    "oat",
]

# This is easier to do directly from the netCDF
fleet_file_pattern = "/disk1/Data/ADS-B/COBALT/20250225_ADS-B.nc"
data = {}
with Dataset(fleet_file_pattern) as ncdf:
    for name in variables:
        data[name] = ncdf.variables[name][:]

# Get wind speeds near Chilbolton at different pressure levels
alt_bins = [
    100 * arguslib.misc.thermo.hPa_to_hft(a) for a in [175, 200, 225, 250, 265, 300]
]

data["dist"] = arguslib.misc.geo.haversine(-1.44, 51.16, data["lon"], data["lat"])
data["alt_bin"] = np.digitize(data["alt_geom"], alt_bins)

output = []
for i in range(24):
    adata = [i]
    for bin in range(7):
        dslice = np.s_[:, (60 * 4 * i) : (60 * 4 * i + 1)]
        mask = (data["dist"][dslice] < 150) & (data["alt_bin"][dslice] == bin)
        if mask.sum() > 0:
            adata.append(np.nanmean(data["ws"][dslice][mask]))
        else:
            adata.append(np.nan)
    output.append(adata)

# Based on comparisons to the Chilbolton wind from OpenMeteo at 1200
# on the 25th Feb, it seems like that this data is indeed in knots.

# The expected windspeed is around 60kph, which is approximately the
# 30 kts seen in the aircraft data. If the wind speeds are too fast in
# the COBALT data, something else is going on here.

# %%
