function h = plottangent(obj, ind, varargin) 
%PLOTTANGENT	Plot waypoints and tangents.
%	PLOTTANGENT(OBJ,IND) plots OBJ.y over OBJ.x, the tangents of
%	the elements of indices IND and highlights the elements IND by
%	a marker '*'.
%	
%	PLOTTANGENT(...,NAME,VALUE) specifies line properties using one
%	or more Name,Value pair arguments.
%	
%	H = PLOTTANGENT(...) returns a handle array H to lineseries
%	objects. H(1,1) is the waypoints-handle, H(i+1,1) and
%	H(i+1,2) the marker- and the tangent-handle of IND(i)
%	respectively.
%	
%	See also WAYPOINTS/PLOT, WAYPOINTS/PLOTDIFF.

	%%% Handle input arguments
	% Check class input ind
	if ~isnumeric(ind)
		error('Input argument IND has to be numeric.');
	end%if

	% Check dimension of input ind
	if (~isvector(ind) && ~isempty(ind))
		error('Input argument IND has to be of size 0xN or Nx0.');
	end%if

	% Apply plot options if unspecified
	if nargin < 3
		opts = {'Color','k'};
	else
		opts = varargin;
	end%if
	opts_obj = {':', 'LineWidth',1};
	opts_marker = {'Marker','x', 'MarkerSize',8};
	opts_tangent = {'--', 'LineWidth',0.5};

	% Initialize handle to graphics objects
	h = gobjects(1+numel(ind), 2);

	% Plot the waypoints and get corresponding axis limtis
	h(1,1) = plot(obj, opts_obj{:}, opts{:});
	xLimits = xlim;
	yLimits = ylim;

	% Check if some tangent indexes are out of range
	if any(ind > numwp(obj))
		warning('WAYPOINTS:plottangent:indxOutOfRange', ...
			['Some indexes exceed the number N=%i of waypoints ', ...
			'These tangents will not be plotted.'], numwp(obj))
		ind(ind > numwp(obj)) = [];
	end%if

	% Plot the tangents and the corresponding waypoints
	hold on
	for i = 1:length(ind)
		iind = ind(i);

		% Marker of tangent point
		h(i+1,1) = plot(obj.x(iind), obj.y(iind), opts_marker{:}, opts{:});

		% Length of tangent
		[r1,r2] = scaleTangentToAxis(xLimits, yLimits, ...
			[obj.x(iind) obj.y(iind)], obj.Head(iind));

		% Start/end point of tangent
		Pstart = [...
			obj.x(iind)+r2*cos(obj.Head(iind)); ...
			obj.y(iind)+r2*sin(obj.Head(iind))];
		Pstop = [...
			obj.x(iind)+r1*cos(obj.Head(iind)); ...
			obj.y(iind)+r1*sin(obj.Head(iind))];

		% Draw the tangent using N arrows sticked together
		N = 25;
		xq = linspace(Pstart(1), Pstop(1), N);
		yq = linspace(Pstart(2), Pstop(2), N);
		% quiver would draw one arrow at every point of xq/yq of
		% length uq/vq; since the last element of xq/yq is the end
		% point, no arrow has to be drawn there
		xq = xq(1:end-1);
		yq = yq(1:end-1);
		uq = ones(size(xq))*(Pstop(1)-Pstart(1))/N;
		vq = ones(size(yq))*(Pstop(2)-Pstart(2))/N;
		scale = 0;
		h(i+1,2) = quiver(xq, yq, uq, vq, scale, opts_tangent{:}, opts{:});

		% Disable legend entries for tangents
		set(get(get(h(i+1, 1), 'Annotation'), 'LegendInformation'), ...
			'IconDisplayStyle','off');
		set(get(get(h(i+1, 2), 'Annotation'), 'LegendInformation'), ...
			'IconDisplayStyle','off');

	end%for
	hold off

	% Set the axis limtis corresponding to waypointss
	axis([xLimits, yLimits]);

end%fcn
