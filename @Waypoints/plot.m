function [h,ax] = plot(varargin)
%PLOT	Plot waypoints.
%	PLOT(OBJ) plots OBJ.y over OBJ.x.
%	
%	PLOT(OBJ,S) additionally applies the line specification S.
%
%	PLOT(...,NAME,VALUE) specifies line properties using one or
%	more Name,Value pair arguments.
%
%	PLOT(AX,...) plots into the axes with handle AX.
%	
%	[H,AX] = PLOT(...) returns the handle H to lineseries objects
%	and axes handle AX.
%	
%	The line specification S is a character string supported by the
%	standard PLOT command. For example
%		PLOT(OBJ, 'LineWidth',2, 'Color',[.6 0 0]) 
%	will plot a dark red line  using a line width of 2 points.
%	
%	NOTE: This method supports non-scalar inputs OBJ!
% 
%	See also PLOT, WAYPOINTS/PLOTG2, WAYPOINTS/PLOTTANGENT,
%	WAYPOINTS/PLOTDIFF.

	% Default plot style if unspecified
	opts = {'LineWidth',1};

	% Handle accepted input argument configurations
	if nargin == 1
		% Easy case: plot() was called with just a WAYPOINTS object
		obj = varargin{1};
		[h,ax] = plot_raw([], obj, opts{:});

	elseif nargin > 1
		% The first argument is either an axes handle or a WAYPOINTS
		% object:
		if isa(varargin{1}, 'matlab.graphics.axis.Axes')
			ax = varargin{1};
			obj = varargin{2};
			if nargin > 2
				opts = varargin(3:end);
			end%if
		else 
			ax = [];
			obj = varargin{1};
			opts = varargin(2:end);
		end%if
		[h,ax] = plot_raw(ax, obj, opts{:});

	end%if

	% Apply plot styles
	grid(ax, 'on');
	axis(ax, 'equal');
	title(ax, getTitleCellString(obj));
	ylabel(ax, 'y [m]');
	xlabel(ax, 'x [m]');

end%fcn
