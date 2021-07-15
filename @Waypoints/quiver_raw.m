function [h,axh] = quiver_raw(axh, obj, varargin)
%QUIVER_RAW		Perform basic quiver-plot operations of waypoints.
%	To be used by public plot methods to avoid multiple calls to
%	plot styles like AXIS, TITLE, XLABEL, ...!
%	
%	NOTE: This method supports non-scalar inputs OBJ!
%
%	See also WAYPOINTS/QUIVER.

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

		u = diff(obj(i).x);
		v = diff(obj(i).y);				
		if strfind([varargin{:}], 'DisplayName')
			h(i) = quiver(axh, obj(i).x(1:end-1), obj(i).y(1:end-1), ...
				u, v, varargin{:}, ...
				'AutoScale','off', ...
				'ShowArrowHead','on');
		else
			name = obj(i).Name;
			if isempty(name)
				name = sprintf('Waypoints %u', i);
			end%if
			h(i) = quiver(axh, obj(i).x(1:end-1), obj(i).y(1:end-1), ...
				u, v, varargin{:}, ...
				'AutoScale','off', ...
				'ShowArrowHead','on', ...
				'DisplayName',name);
		end%if
	end%for

	if numel(obj) > 1
		legend('show', 'Location','best');
	end%if

	% Reset axes to initial state
	set(axh, 'NextPlot',npState);

end%fcn
