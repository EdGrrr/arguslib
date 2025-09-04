.. _user_guide:

User Guide
==========

Welcome to the ``arguslib`` user guide. 
This guide provides an overview of the core concepts and classes for visualizing and analyzing coincident data from cameras, radar, and other instruments.

The Core Abstractions
---------------------

Arguslib is made up of geolocation and visualization components, all based around structures representing instruments.
The two most fundamental concepts are the ``Position``, which defines a location in space, and the ``PlottableInstrument``, which defines an object that can be visualized.

The `Position` Class: The Foundation of Geolocation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fundamental concept for geolocation is the ``Position`` class.
It stores the longitude, latitude, and altitude of a point and provides methods for calculating distances and bearings to other points. 
It is the fundamental building block for locating physical instruments and annotating data.

.. autoclass:: arguslib.instruments.instruments.Position
   :members: target_ead, target_xyz, ead_to_lla

The `PlottableInstrument`: A Universal Drawing API
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once you have a location, you need something to visualize it on. The ``PlottableInstrument`` is an abstract class that defines a standard interface for any object that can be visualized. It guarantees that any ``PlottableInstrument`` will have ``show()`` and ``annotate_positions()`` methods.

This is a powerful concept because it applies not only to physical devices, but also to higher-level "interfaces" that combine plottables and other data sources.

.. autoclass:: arguslib.instruments.instruments.PlottableInstrument
   :members: show, annotate_positions

The `Instrument`: Physical Devices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``Instrument`` class represents a physical device, like a camera or a radar. It extends ``PlottableInstrument`` by adding the concept of a physical location (``position``) and orientation (``rotation``) in the world.
It also introduces new instrument-relative coordinate systems including the ``iead`` (instrument elevation, azimuth and distance, measured relative to some reference axes, such as the image plane and focal axes).

.. autoclass:: arguslib.instruments.instruments.Instrument
   :members:
   :undoc-members:
   :exclude-members: show, annotate_positions, get_data_time, initialise_data_loader

Data Loaders
~~~~~~~~~~~~

Each ``Instrument`` has a ``data_loader`` attribute. This object is responsible for finding and loading the actual data (e.g., image frames or radar scans) for a given time. This separation of concerns keeps the instrument logic clean from the details of file formats and directory structures. The data loader is typically initialized automatically when data is first requested.

Concrete Instruments: The Camera
--------------------------------

The ``Camera`` class is a concrete implementation of an ``Instrument``. It represents a single camera, handles its specific calibration, and knows how to load and display image frames.

.. autoclass:: arguslib.camera.camera.Camera
   :members: from_config, show, annotate_positions
   :show-inheritance:

Specialized Cameras
~~~~~~~~~~~~~~~~~~~

The power of this structure is that it's easy to extend. For example, ``arguslib`` provides specialized camera classes that build upon the base ``Camera``:

*   **UndistortedCamera**: Inherits from ``Camera`` but overrides the data loading method to apply a fisheye undistortion to the image before displaying it.
*   **DirectCamera**: A more advanced subclass that bypasses Matplotlib entirely, drawing annotations directly onto the image array using OpenCV. This is ideal for generating videos efficiently.

.. autoclass:: arguslib.camera.undistorted_camera.UndistortedCamera
   :members:
   :show-inheritance:

.. autoclass:: arguslib.camera.direct_camera.DirectCamera
   :members: show, annotate_positions, image
   :show-inheritance:

Combining Instruments and Data: Interfaces
------------------------------------------

Interfaces are also ``PlottableInstrument`` objects, but they are not physical devices. Instead, they **compose** one or more other plottable instruments to create more complex, synchronized visualizations. This is a key design pattern in ``arguslib``.

The ``RadarInterface``
~~~~~~~~~~~~~~~~~~~~~~

The ``RadarInterface`` takes a ``Radar`` and another ``PlottableInstrument`` (like a ``Camera`` or even a ``CameraArray``) and displays them side-by-side.

.. autoclass:: arguslib.radar.radar_interface.RadarInterface
   :members: from_campaign, show, annotate_positions
   :undoc-members:
   :show-inheritance:

The ``AircraftInterface``
~~~~~~~~~~~~~~~~~~~~~~~~~

Similarly, the ``AircraftInterface`` wraps another ``PlottableInstrument`` and adds the capability to load and display aircraft flight tracks on top of it. Because it is itself a ``PlottableInstrument``, it can be passed to other interfaces for even more complex plots.

.. autoclass:: arguslib.aircraft.aircraft_interface.AircraftInterface
   :members: from_campaign, show, plot_trails
   :show-inheritance:
