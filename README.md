# Waypoints

A MATLAB class to work with data representing waypoints in a cartesian coordinate frame.

## Introduction
To create an instance of the `Waypoints` class use one of the following approaches:
- Call the class constructor `Waypoints`.
- Call static method `xy2Waypoints` providing x/y coordinates.
- Call static method `ll2Waypoints` providing longitudinal/lateral coordinates, e.g. from GPS.
- Call static method `pp2Waypoints` providing a piecewise polynomial structure.

## Methods
The `Waypoints` class provides multiple methods for modification and visualization.
For example shift, rotate, append or reorder waypoints.
For a list of available methods type `help Waypoints` or `methods Waypoints` in your MATLAB command window.
