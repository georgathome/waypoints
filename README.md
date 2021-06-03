# waypoints

A MATLAB class to work with data representing cartesian waypoints.

## Introduction
To create an instance of the `Waypoints` class use one of the folling approaches:
- Call the class constructor `Waypoints`.
- Call static method `xy2Waypoints` providing x/y coordinates.
- Call static method `ll2Waypoints` providing longitudinal/lateral GPS coordinates.
- Call static method `pp2Waypoints` providing a piecewise polynomial structure.

## Methods
The `Waypoints` class provides multiple methods for modification and visualization.
For example shift, rotate, append or reorder waypoints.
For a list of available methods type `help Waypoints` or `methods Waypoints` in your MATLAB command window.
