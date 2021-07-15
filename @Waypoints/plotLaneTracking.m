function plotLaneTracking(obj_global,...
		xyCG_global, psi_global, LAD, indx, sd_T, xyCG_T)


	%%% handle input arguments
	% number of input arguments
	if nargin ~= 5 && nargin ~= 7
		error('Invalid number of input arguments!');
	end%if

	% XYCG_GLOBAL
	xyCG_global = xyCG_global(:);
	if any(size(xyCG_global) ~= [2 1])
		error('Input argument XYCG_GLOBAL must be a vector of two elements!');
	end%if

	% PSI_GLOBAL
	if ~isscalar(psi_global)
		error('Input argument PSI_GLOBAL must be a scalar!');
	end%if

	% LAD
	LAD = LAD(:);
	if ~isvector(LAD)
		error('Input argument LAD must be a vector!');
	end%if


	%%% start plotting
	% desired waypoints (global)
	h = plottangent(obj_global,indx);
	set(h(1),...
		'LineStyle','none',...
		'Marker','o',...
		'MarkerSize',3,...
		'MarkerFaceColor','b');

	% vehicle representation (global)
	hold on
	plot(xyCG_global(1),xyCG_global(2),...
		'ob',...
		'MarkerSize',7,...
		'MarkerFaceColor','b');
	Fzglachse = xyCG_global + max(LAD)*[cos(psi_global); sin(psi_global)];
	plot([xyCG_global(1),Fzglachse(1)],[xyCG_global(2),Fzglachse(2)],...
		'b',...
		'LineWidth',2)


	if nargin > 5

		%%% handle input arguments
		% SD_T
		if ~isa(sd_T, 'Waypoints')
			warning('WAYPOINTS:plotLaneTracking:class',...
				'Input argument SD_T not of class WAYPOINTS!');
% 					return;
			disp('Try to convert to class WAYPOINTS...')
			sd_T = Waypoints(...
				sd_T.x,...
				sd_T.y,...
				sd_T.s,...
				sd_T.Curv,...
				sd_T.Head,...
				sd_T.Type,...
				ones(size(sd_T.x)));
			disp('... done!');
		end%if

		% XYCG_T
		xyCG_T = xyCG_T(:);
		if any(size(xyCG_T) ~= [2 1])
			error('Input argument XYCG_T must be a vector of two elements!');
		end%if


		%%% start plotting
		% Solltrajektorie (transformiert)
		h_V = plottangent(sd_T,indx, 'Color','k');
		set(h_V(1),...
			'Color','g',...
			'LineStyle','none',...
			'Marker','x',...
			'MarkerSize',5,...
			'MarkerFaceColor','r');

		% Fahrzeug (transformiert)
		hold on
		plot(xyCG_T(1),xyCG_T(2),...
			'rx','MarkerSize',7,...
			'MarkerFaceColor','r');
		Fzglachse_T = xyCG_T + [max(LAD);0];
		plot([xyCG_T(1),Fzglachse_T(1)],[xyCG_T(2),Fzglachse_T(2)],...
			'r',...
			'LineWidth',2)

	end%if

	hold off
	grid on
	axis auto
	axis equal

end%fcn
