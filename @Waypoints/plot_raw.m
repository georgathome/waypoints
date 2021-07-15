function [h,axh] = plot_raw(axh, obj, varargin)
%PLOT_RAW	Perform basic plot operations of waypoints.
%	To be used by public plot methods to avoid multiple calls to
%	plot styles like AXIS, TITLE, XLABEL, ...!
%
%	NOTE: This method supports non-scalar inputs OBJ!
% 
%	See also WAYPOINTS/PLOT.

	% Get current status of axes 'NextPlot' property
	if isempty(axh) || ~ishghandle(axh)
		axh = gca;
	end%if
	npState = get(axh, 'NextPlot');

	% Plot paths
	h = gobjects(numel(obj), 1);
	for i = 1:numel(obj)
		if i == 2
			set(axh, 'NextPlot', 'add');
		end%if


		if strfind([varargin{:}], 'DisplayName')
			h(i) = plot(axh, obj(i).x, obj(i).y, varargin{:});
		else
			name = obj(i).Name;
			if isempty(name)
				name = sprintf('Waypoints %u', i);
			end%if
			h(i) = plot(axh, obj(i).x, obj(i).y, varargin{:}, ...
				'DisplayName',name);
		end%if
	end%for

	if numel(obj) > 1
		legend(axh, 'show', 'Location','best');
	end%if

	% Reset axes to initial state
	set(axh, 'NextPlot',npState);

end%fcn
