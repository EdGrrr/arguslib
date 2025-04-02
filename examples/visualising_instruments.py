# %%
from matplotlib import pyplot as plt
from arguslib.instruments.camera import Camera
import datetime as dt

from arguslib.instruments.radar import Radar

# %% Test default camera logic (on polar plot, with the axes aligned with the ordinal directions)
cam = Camera.from_config("COBALT", "3-8")
cam.show(dt.datetime(2025, 3, 25, 15))


# %% We can also plot the image "unflipped", either by setting up the axes to be flipped:
cam.show(dt.datetime(2025, 3, 25, 15), theta_behaviour="unflipped_ordinal_aligned")

# %% Or by setting the lr_flip argument to False
cam.show(
    dt.datetime(2025, 3, 25, 15),
    lr_flip=False,
)

# %% We can set the theta behaviour to be "pixels", here the theta grid shows the major axes of the pixel grid.
cam.show(dt.datetime(2025, 3, 25, 15), theta_behaviour="pixels")

# %% This lines up with the native image grid, which is seen when we plot on non-polar axes
fig, ax = plt.subplots()
cam.show(dt.datetime(2025, 3, 25, 15), ax=ax)

# %% To still get the rotation, but in place of an existing non-polar set of axes, we can use the replace_ax argument
fig, ax = plt.subplots()
cam.show(dt.datetime(2025, 3, 25, 15), replace_ax=ax)

# %% The infrastructure can be rigged to plot in a way that doesn't make any sense
# here, we set up a "bearing" axes, but then don't flip the axes, so the angles don't line up with the ordinal directions.
cam.show(
    dt.datetime(2025, 3, 25, 15),
    theta_behaviour="bearing",
    lr_flip=False,
)

# %% Radar objects can also be plotted using the analogous show function.
# The axes don't have so much rotating and flipping to do.
# Following pyart default behaviour, the positive x direction is towards 0 degrees azimuth. (i.e. south to north).
radar = Radar.from_config("COBALT")
radar.show(dt.datetime(2025, 3, 25, 9), var="DBZ")


# %%
