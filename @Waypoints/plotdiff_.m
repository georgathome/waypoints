function h = plotdiff_(obj, fh, prop)
%PLOTDIFF_	Plot WAYPOINTS object property with specific appearance.
%	PLOTDIFF_(OBJ,FH,PROP) works similar to PLOTDIFF but plots the
%	property PROP.
%	
%	PROP can take the following values:
%	 'curv': plots the curvature over normed index
%	 'head': plots the angle over normed index
%	 's':	 plots curve length over normed index
%	 'type': plots waypoints type over normed index
%	 'nbr':	 plots segment number over normed index
%	
%	FH can be specified as []. In this case, the default value is
%	used.
%	
%	See also WAYPOINTS/PLOTDIFF.

	N = numwp(obj);
	xx = (0:N-1)/N;
	xLblString = 'index [1/index_{max}]';

	switch lower(prop)
		case 'curv'
			yy = obj.Curv;
			yLblString = 'curvature [1/m]';

		case 'head'
			yy = obj.Head;
			yLblString = 'heading [rad]';

		case 's'
			yy = obj.s;
			yLblString = 's [m]';

		case 'type'
			yy = obj.Type;
			yLblString = 'type [-]';

		case 'nbr'
			yy = obj.Nbr;
			yLblString = 'nbr [-]';

		otherwise
			error(['Unknown property string. ',...
				'Type help plotdiff_ for valid property strings.']);

	end%switch

	% create WAYPOINTS object for plotting
	obj_new = Waypoints(xx, yy, obj.s, obj.Curv, obj.Head,	obj.Type, obj.Nbr);

	% plot
	h = plotdiff(obj_new, fh);

	% apply plot styles
	axis normal
	grid on
	xlabel(xLblString);
	ylabel(yLblString);
	title('');

end%fcn
