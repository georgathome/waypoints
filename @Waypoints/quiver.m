function h = quiver(obj, varargin)
%QUIVER		Quiver plot of waypoints.
%	QUIVER(OBJ) plots vectors as arrows at coordinates
%	(OBJ.y,OBJ.x) with components (diff(OBJ.x),diff(obj.y))
%	
%	QUIVER(OBJ,S) additionally applies the line specification
%	S. Property 'AutoScale' is disabled permanently!
%	
%	H = QUIVER(...) returns the handle H to Quiver object.
%	
%	The line specification S is a character string supported by the
%	standard QUIVER command.
%	
%	See also QUIVER.


	% apply plot options if unspecified
	if isempty(varargin)
%				varargin = {'o','MarkerSize',2,'MarkerFaceColor','blue'};
		varargin = {'-b', 'LineWidth',1};
	end%if

	% Plot waypoints
	h = quiver_raw([], obj, varargin{:});

	% apply plot styles
	grid on;
	axis equal;
	title(getTitleCellString(obj));
	ylabel('y [m]');
	xlabel('x [m]');

end%fcn
