% Testfile: 
% 
% Tests for ...

clear;
close all
clc;


%% Load data
gps	= importdata('gpsrec.mat');


%% Plot
p = geoplot(gps.GPS_Latitude.Data, gps.GPS_Longitude.Data,'ro'); 
geobasemap('satellite')


%% Resample data set
% Create PATH2D object from Lat/Lon
r = Waypoints.ll2Waypoints(gps.GPS_Latitude.Data, gps.GPS_Longitude.Data, ...
	'head', gps.GPS_Heading.Data*pi/180);
r = shiftTo(r);
r.Name = 'RAW';

% Resample using different methods
deltaSet = 2;
rr1 = resample(r, deltaSet, 'pchip');
rr1.Name = 'pchip';
rr2 = resample(r, deltaSet, 'spline');
rr2.Name = 'spline';

% Visualize results
hold off
[~, ax] = curvplot(r, 'ro');
set(ax, 'NextPlot', 'add');
curvplot(ax, rr1,'b.');
curvplot(ax, rr2,'g.');

% Spline interpolation/approximation
ppxy = spline(r.s, [r.x, r.y]');
rs = Waypoints.pp2Waypoints(0:1:r.s(end), ppxy);
rs.Name = 'interp.';
curvplot(ax, rs, 'k--');
set(ax, 'Nextplot','replace');

legend(ax(1), 'show')
