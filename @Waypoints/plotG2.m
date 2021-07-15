function [h,ax] = plotG2(varargin)
%PLOTG2		Plot waypoints, heading and curvature.
%
%	For Syntax see:
%	See also WAYPOINTS/PLOT.

	% Handle accepted input argument configurations
	if nargin == 1
		% Easy case: plot() was called with just a WAYPOINTS object 
		ax = [];
		obj = varargin{1};
		opts = {};

	elseif nargin > 1
		% The first argument is either an axes handle or a WAYPOINTS
		% object:
		if isa(varargin{1}, 'matlab.graphics.axis.Axes')
			ax = varargin{1};
			obj = varargin{2};
			opts = varargin(3:end);
		else 
			ax = [];
			obj = varargin{1};
			opts = varargin(2:end);
		end%if

	end%if

	if isempty(ax) || any(~ishghandle(ax))
		ax = [subplot(2, 2, [1,3]); subplot(2, 2, 2); subplot(2, 2, 4)];
	elseif numel(ax) < 3
		error('If specified, three axes handles are required!');
	end%if

	N = numel(obj);
	h = gobjects(N,3);
	h(:,1) = plot(ax(1), obj, opts{:});
% 			npStatus = get(ax(2:3), 'NextPlot');
% 			set(ax(2:3), 'NextPlot','replace');
	for i = 1:N
		if i == 2
			set(ax(2:3), 'NextPlot','add');
		end
		h(i,2) = plot(ax(2), obj(i).s, obj(i).Head, opts{:});
		h(i,3) = plot(ax(3), obj(i).s, obj(i).Curv, opts{:});
	end%for
% 			set(ax(2:3), {'NextPlot'},npStatus);

	grid(ax(2), 'on');
	ylabel(ax(2), 'Heading [rad]');
	grid(ax(3), 'on');
	xlabel(ax(3), 'Length [m]');
	ylabel(ax(3), 'Curvature [1/m]');

	% linkaxes(_, 'x') sets axes XLimMode to manual which can cause
	% data of subsequent calls of this method to be not completely
	% shown. Therefore reset axes limit mode to 'auto'.
	linkaxes(ax(2:3), 'x'); 
	set(ax, 'XlimMode', 'auto');

end%fcn
