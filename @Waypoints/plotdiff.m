function h = plotdiff(obj, fh)
%PLOTDIFF	Plot waypoints with specific appearance.
%	PLOTDIFF(OBJ) plots the different types of waypoints OBJ using
%	the according pre-defined type-specific color and marker.
%	
%	PLOTDIFF(OBJ,FH) marker symbols are switched periodically at
%	indices specified by the function handle FH. Some usefull
%	values of FH might be
%	  (D) @(OBJ) diff(OBJ.Type) 
%	  (.) @(OBJ) diff(OBJ.Type) | diff(sign(OBJ.Curv))
%	  (.) @(OBJ) diff(OBJ.Nbr)
%	where (D) is used by PLOTDIFF(OBJ). You can also
%	specify the indices where marker symbols change manually by
%	using something like
%		.) @(OBJ) [0 0 1 0 1 0 0 0 1 ...].
%	Internally FH is evaluated to FH(OBJ) ~= 0.
%	
%	H = PLOTDIFF(...) returns a column vector of handles to
%	lineseries objects, one handle per WAYPOINTS type.
%	
%	See also WAYPOINTS/PLOT, WAYPOINTS/PLOTTANGENT.


	% handle input arguments
	if nargin < 2 || isempty(fh)
		fh = @(arg) diff(arg.Type); 
	end%if

	if ~isa(fh, 'function_handle')
		error('Second input argument must be of class function handle!')
	end%if


	% index numbering of number of elements of OBJ
	indBase = 1:numwp(obj);


	% get the segment-grouping indices and the number of indices
	flag = fh(obj) ~= 0;
	if length(flag) > length(indBase)
		error(['The length %u of your evaluated function handle ',...
			'FH exceeds the number %u of waypoints.'],...
			length(flag), length(indBase));
	end%if
	ind = [0, indBase(flag), indBase(end)];
	nbrInd = length(ind);


	% define/extend the plot settings to the number of segments to
	% plot
	colorStyles = {'m','b','g','r'};
	markerStyles = {...
		'o','diamond';... % Plot marker marker symbol
		5,	4}; % Plot marker marker size
	n = ceil(nbrInd/size(markerStyles, 2));
	plotMarker_ = repmat(markerStyles, 1, n);


	% get axis current 'hold on' status which might be set from
	% outside
	axh = gca;
	npState = get(axh, 'NextPlot');

	%%% plot the segments with according plot options
	h = gobjects(nbrInd-1, 1);
	for i = 1:nbrInd-1

		% get the segment data to plot
		sd = getSamples(obj, ind(i)+1:ind(i+1));

		if i == 2
			hold on
		end
		h(i) = plot_raw([], sd,... % TODO: pass axes handle?
			'LineStyle','none',...
			'Color',	colorStyles{sd.Type(1)+2},...
			'Marker',	plotMarker_{1,i},...
			'MarkerSize', plotMarker_{2,i});

	end%for

	% reset axis 'hold on' status (might have been set from outside)
	set(axh, 'NextPlot',npState);

	% unsure about usefulness
	% line style plotting if no additional ....
% 			if nargin < 2
% 				set(h,'LineStyle','-','LineWidth',2,'Marker','none');
% 			end%if

	% apply plot styles
	grid on;
	axis equal;
	title(getTitleCellString(obj));
	ylabel('y [m]');
	xlabel('x [m]');

end%fcn
