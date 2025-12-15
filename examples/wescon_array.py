import datetime
from arguslib.camera.camera_array import CameraArray
from arguslib.instruments.instruments import Position
from arguslib.misc.plotting import show_the_sun
import matplotlib.pyplot as plt

dt = datetime.datetime(2023, 8, 30, 9, 30, 0)
basic_array = CameraArray.from_config("WesConArray")
axes = basic_array.show(dt)

# Plot a column of points above Camra/Chilbolton big dish
basic_array.annotate_positions(
    [Position(-1.4419, 51.1553, i) for i in range (10)],
    dt,
    ax=axes,
    marker="o",
    lw=0,
    color="darkred",
    label="Chilbolton",
)
plt.show()
