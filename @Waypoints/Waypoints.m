classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) Waypoints
%WAYPOINTS	Waypoints in 2-D cartesian space.
%	
%	OBJ = WAYPOINTS(X,Y,S,K,PHI,TYPE,NBR) creates an object OBJ of class
%	WAYPOINTS representing 2-D waypoints using the properties listed below.
% 
%	OBJ = WAYPOINTS() creates an instance OBJ of class WAYPOINTS using
%	default values.
%	
%	A curvature K > 0 (K < 0) indicates a left (right) turn by moving from
%	[X(1),Y(1)] -> [X(end),Y(end)].
%	
%	NBR is usually 1, but becomes interesting when multiple WAYPOINTS are
%	connected together. Then NBR represents the number of the WAYPOINTS
%	object in the list of connected WAYPOINTS.
%	
%	WAYPOINTS Properties:
%	 NAME	- Name of the waypoints.
%	 X		- x-coordinate [m].
%	 Y		- y-coordinate [m].
%	 S		- Length parameter, strictly increasing.
%	 HEAD	- Heading angle [rad].
%	 CURV	- Curvature [1/m].
%	 TYPE	- Integer indicating the waypoint type.
%	 NBR	- Segment number of joined WAYPOINTS.
%	 INITPOINT - Initial point [x(1) y(1)].
%	 TERMPOINT - Terminal Point [x(end) y(end)].
%	
%	WAYPOINTS Methods:
%	 - INSTANTIATION
%	   WAYPOINTS - Class constructor.
%	   ll2Waypoints - Create WAYPOINTS from latitude/longitude coordinates.
%	   pp2Waypoints - Create WAYPOINTS from piecewise polynomial structure.
%	   xy2Waypoints - Create WAYPOINTS from x/y coordinates.
%	
%	 - MANIPULATION
%	   append	- Append waypoints objects.
%  	   changeSignOfCurvature - Change sign of curvature property CURV.
%	   getSamples - Select subset of waypoints.
%	   perpendicularDistance - Perpendicular distance from line to waypoints.
%	   plus		- Connect waypoints using '+'.
%	   resample - Resample waypoints to finer/coarser resolution.
%	   reverse	- Reverse waypoints ordering.
%	   rotate	- Rotate waypoints by an angle.
%	   setStartIndex - Set waypoint to be the initial waypoint.
%	   shiftBy	- Shift waypoints by offset.
%	   shiftTo	- Shift waypoints to point.
%	   unwrap	- Unwrap property HEAD.
% 
%	 - POINT REDUCTION/CURVE FITTING
%	   numwp		 - Number of waypoints.
%	   douglasPeuker - Ramer�Douglas�Peucker point reduction algorithm. 
%	   fitStraight	 - Fit straight line to waypoints.
%	   fitCircle	 - Fit circle to waypoints.
%	 
%	 - VISUALIZATION
%	   plot		 - Plot waypoints.
%	   plotdiff	 - Plot waypoints with specific appearance.
%	   plotdiff_ - Plot waypoints properties with specific appearance.
%	   plotG2	 - Plot waypoints, heading and curvature.
%	   plottangent - Plot waypoints and tangents.
%	   quiver	 - ToDo
%	 
%	 - EXPORT
%	   write2file - Write WAYPOINTS properties to file.
% 
%	 - MISC
%	   laneTracking - Get the lane tracking pose.
%

% DEVELOPMENT NOTES:
% 
%	(1) Add methods for plot labeling (x/y labels, title, ...). DONE
% 
%	(2) Show the clothoid direction in plots.
% 
%	(3) For implementing methods that accept a graphics object its first
%	argument (for example, an axes handle):
%	https://www.mathworks.com/help/matlab/matlab_oop/overloading-plotting-functions.html

% Subject: 2-D Waypoints.
% Author: $Author$
% Date: $LastChangedDate$
% Revision: $Revision$

	
	properties
		% NAME - The name of the WAYPOINTS object (char or string).
		Name = '';
	end%properties
	
	properties (SetAccess = private)
		% X - X-coordinates [m].
		%	n-by-1
		x = [0; 1];
		
		% Y - Y-coordinates [m].
		%	n-by-1
		y = [0; 0];
		
		% S - Cumulative length of consecutive waypoints [m].
		%	n-by-1
		s = [0; 1];
		
		% HEAD - Waypoint heading angle [rad].
		%	n-by-1
		%	Measured from the positive x-axis to the waypoint's tangent.
		Head = [0; 0];
		
		% CURV - Waypoint curvature [1/m].
		%	n-by-1
		Curv = [0; 0];
		
		% TYPE - User defined waypoint type.
		%	n-by-1 int8
		Type = zeros(2, 1, 'int8');
		
		% NBR - Number of waypoint segment.
		%	n-by-1 uint16
		%	When appending multiple WAYPOINTS objects to one single
		%	WAYPOINTS object, this property allows to distinguish them in
		%	the single object.
		%
		%	Needs to be monotonically increasing and start with 1!
		Nbr = ones(2, 1, 'uint16');
	end%properties
	
	properties (Dependent)
		% INITPOINT - The initial waypoint.
		InitPoint
		
		% TERMPOINT - The terminal (end) waypoint.
		TermPoint
	end%properties
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	methods
		
		function obj = Waypoints(x, y, s, curv, head, type, nbr) %#codegen
		%WAYPOINTS	Create an instance of the WAYPOINTS class.
			
			% Allow class instantiation with no input arguments by using
			% default values
			if nargin < 1
				return
			end%if
			
			%%% (1) Check input argument sizes
			% Properties X, Y, S, HEAD and CURV are required to have the
			% same size
			Nx = numel(x);
			assert(isequal(Nx, numel(y), numel(s), numel(head), numel(curv)), ...
				'WAYPOINTS:Waypoints:numelInputs',...
				'Input arrays must have the same number of elements!');
			
			% Before we can check the size of optional input arguments TYPE
			% and NBR, we must ensure they are defined
			if nargin < 7
				Nbr = ones(Nx, 1, 'uint16');
			elseif isscalar(nbr)
				% Explicit array size syntax for compatibility in MATLAB
				% Function blocks within Simulink.
				Nbr = nbr(1) * ones(Nx, 1, class(nbr));
			else
				Nbr = nbr;
			end%if
			if nargin < 6
				Type = -ones(Nx, 1, 'int8');
			elseif isscalar(type)
				Type = type(1) * ones(Nx, 1, class(type));
			else
				Type = type;
			end%if
			
			assert(isequal(numel(x), numel(Nbr), numel(Type)), ...
				'WAYPOINTS:Waypoints:numelTypeNbrInputs',...
				'Size mismatch for input TYPE and/or NBR!');
			
			
			%%% (2) Check for NaN's in waypoints
			isNotNan = ~or(isnan(x), isnan(y));
			if any(~isNotNan)
				fprintf('Warning: NaN''s samples have been removed!\n');
			end%if
			
			
			%%% (3) Set class properties
			obj.x		= reshape(x(isNotNan), [], 1);
			obj.y		= reshape(y(isNotNan), [], 1);
			obj.s		= reshape(s(isNotNan), [], 1);
			obj.Head	= reshape(head(isNotNan), [], 1);
			obj.Curv	= reshape(curv(isNotNan), [], 1);
			obj.Type	= reshape(Type(isNotNan), [], 1);
			obj.Nbr		= reshape(Nbr(isNotNan), [], 1);
			
		end%CONSTRUCTOR
        
		function obj = append(obj, varargin)
		%APPEND		Concatenate waypoints.
		%	OBJ = APPEND(OBJ1,OBJ2,...,OBJN) appends waypoints OBJ2 at the
		%	end of waypoints OBJ1 and so on to create waypoints OBJ.
		%	
		%	In contrast to PLUS, APPEND does not shift any of the
		%	waypoints. Also, the total number of waypoints is equal to the
		%	sum of waypoints of OBJ1 to OBJN.
		%
		% See also WAYPOINTS/PLUS.
			
			if nargin < 2
				% return single input argument
				return;
			end%if
			
			% Add WAYPOINTS from VARARGIN iteratively
			for i = 1:numel(varargin)
				obj2append = varargin{i};
				
				% Calculate the distance from OBJ.TERMPOINT to OBJ2.INITPOINT
				ds = sqrt( sum((obj2append.InitPoint - obj.TermPoint).^2) );

				obj = Waypoints(...
					[obj.x;		obj2append.x], ...
					[obj.y;		obj2append.y], ...
					[obj.s;		obj2append.s + obj.s(end) + ds], ...
					[obj.Curv;	obj2append.Curv], ...
					[obj.Head;	obj2append.Head], ...
					[obj.Type;	obj2append.Type], ...
					[obj.Nbr;	obj2append.Nbr + obj.Nbr(end)]);
			end%for
			
		end%fcn
		
		function obj = changeSignOfCurvature(obj)
		%CHANGESIGNOFCURVATURE	Flip waypoints curvature sign.
		%	OBJ = CHANGESIGNOFCURVATURE(OBJ) changes the curvature sign of
		%	waypoints. Properties Y and HEAD are modified accordingly while
		%	maintaining the initial point.
		%	
		%	NOTE: This method supports array inputs OBJ!
		
			for i = 1:numel(obj)
				% Get starting point/angle
				P0 = obj(i).InitPoint;
				phi0 = obj(i).Head(1);
			
				% Shift to origin and rotate so initial slope is zero
				obj(i) = shiftTo(obj(i));
				obj(i) = rotate(obj(i), -phi0);
			
				% Invert the curvature
				obj(i) = Waypoints(...
					+obj(i).x,...
					-obj(i).y,...
					+obj(i).s,...
					-obj(i).Curv,...
					-obj(i).Head,...
					+obj(i).Type,...
					+obj(i).Nbr);
			
				% Undo the shift/rotate procedure
				obj(i) = rotate(obj(i), phi0);
				obj(i) = shiftTo(obj(i), P0);
			end%for
			
		end%fcn
		
		function [obj,idx] = douglasPeuker(obj, eps)
		%DOUGLASPEUKER	Ramer-Douglas-Peucker point reduction.
		%	OBJR = DOUGLASPEUKER(OBJ,EPS) applies the Ramer-Douglas-Peuker
		%	point reduction algorithm with parameter EPS. None of the
		%	removed waypoints has a distance greater EPS to the resulting
		%	waypoints!
		%	
		%	[...,IDX] = DOUGLASPEUKER(OBJ,EPS) returns an array IDX so that
		%	OBJR = GETSAMPLES(OBJ, IDX).
		%
		%	NOTE: This implementation was inspired by dpsimplify.m by
		%	Wolfgang Schwanghart found at MathWorks File Exchange.
		% 
		%	See also WAYPOINTS/GETSAMPLES.
			
			% Initialize a logical array indicating which waypoints to keep
			N = numwp(obj);
			keepIdx = true(1, N);
			
			% Recursively set indexes of waypoints that can be discarded to
			% false
			dprec(1, N);
			
			function dprec(idx0, idx1)
				d = perpendicularDistance(obj.getSamples(idx0:idx1));
				[val_max,idx_max] = max(abs(d));
				if val_max > eps
					% Split waypoints at IDX_SPLIT and call recursion with
					% those two resulting segments until we end in the else
					% statement.
					idx_split = idx_max + idx0 - 1;
					dprec(idx0, idx_split);
					dprec(idx_split, idx1);
				else
					 if idx0 ~= idx1-1
						keepIdx(idx0+1:idx1-1) = false;
					 end%if
				end%if
			end%fcn
			
			idx = find(keepIdx);
			obj = obj.getSamples(idx);
			
		end%fcn
		
		function [straight,e,k,d] = fitStraight(obj, indMinMax, doPlot)
		%FITSTRAIGHT	Fit a straight line to waypoints.
		%	[STRAIGHT,E,K,D] = FITSTRAIGHT(OBJ) fits a straight line with
		%	slope K and offset D to a set of given waypoints (OBJ.X, OBJ.Y)
		%	minimizing the error E. Returns STRAIGHT of class WAYPOINTS .
		%	
		%	[...] = FITSTRAIGHT(OBJ,INDMINMAX) lets you specify a start
		%	index I = INDMINMAX(1) and end index U = INDMINMAX(2) for
		%	considering only the range I:U for the fitting procedure.
		%	
		%	[...] = FITSTRAIGHT(OBJ,INDMINMAX,DOPLOT) allows to disable the
		%	plot for checking the fitting result visually, which is enabled
		%	by default.
		%	
		%	Parameters C and D model the fitted line according to
		%	  y_fit = C*OBJ.X + D
		%	
		%	Minimization is performed in the least-squares sense minimizing
		%	the sum of squared errors:
		%	  SUM[(Y(i)-C*X(i) - D)^2]
		%	   i
		%	
		%	See also WAYPOINTS/FITCIRCLE.
			
			% Handle input arguments
			if (nargin < 2) || isempty(indMinMax)
				indMin = 1;
				indMax = numwp(obj);
			else
				indMin = indMinMax(1);
				indMax = indMinMax(2);
			end%if
			
			if nargin < 3
				doPlot = true;
			end%if
			
			% Extract relevant x/y data
			xsub = obj.x(indMin:indMax);
			ysub = obj.y(indMin:indMax);
			
			% Create (overdetermined) system of equations
			A = [xsub, ones(indMax-indMin+1, 1)];
			b = ysub;
			
			% Solve system of equations: A*[c;d] = b, where y = c*x+d
			cd = (A'*A)\A'*b;
			k = cd(1);
			d = cd(2);
			
			% The fitted y coordinates
			yfit = k*xsub + d;
			
			% Create WAYPOINTS object
			dx = xsub(end) - xsub(1);
			dy = yfit(end) - yfit(1);
			straight = Waypoints(xsub, yfit, ...
				sFrom_x_y(xsub, yfit), ...
				zeros(1, numel(xsub)), ...
				atan2(dy, dx)*ones(numel(xsub), 1), ...
				0, ...
				1);
			
			% Calculate error
			e = 1/numel(ysub)*sum((ysub - yfit).^2);
			
			% Plot if requested
			if doPlot
				plot(obj, 'r', 'DisplayName','Waypoints');
				hold on
				plot(straight, 'b', 'DisplayName','Straight');
				hold off
				legend('show');
			end%if
			
		end%fcn
		
		function [circle,e,xc,yc,R] = fitCircle(obj, N, indMinMax, doPlot)
		%FITCIRCLE	Fit a circle to waypoints.
		%	[CIRCLE,E,XC,YC,R] = FITCIRCLE(OBJ,N) fits a circle with
		%	center(XC,YC) and radius R to waypoints (OBJ.X, OBJ.Y)
		%	minimizing the error E. Returns CIRCLE of clas WAYPOINTS
		%	sampled at angles of linspace(0, 2*pi, N).
		%	
		%	[...] = FITCIRCLE(...,INDMINMAX) lets you specify a start index
		%	I = INDMINMAX(1) and end index U = INDMINMAX(2) for considering
		%	only the range I:U for the fitting procedure.
		%	
		%	[...] = FITCIRCLE(...,INDMINMAX,DOPLOT) allows to disable the
		%	plot for checking the fitting result visually, which is enabled
		%	by default.
		%	
		%	Minimization is performed in the least-squares sense minimizing
		%	the sum of squared errors:
		%	  SUM[(R(i)^2-R^2)^2]
		%	   i
		%	
		%	See also WAYPOINTS/FITSTRAIGHT.
			
			% Handle input arguments
			if nargin < 3 || isempty(indMinMax)
				indMin = 1;
				indMax = numwp(obj);
			else
				indMin = indMinMax(1);
				indMax = indMinMax(2);
			end%if
			
			if indMin >= indMax
				error('WAYPOINTS:FITCIRCLE:index',...
					'Start index must be smaller than end index!');
			end%if
			
			if nargin < 4
				doPlot = true;
			end%if
			
			
			% Extract relevant x/y data
			xsub = obj.x(indMin:indMax);
			ysub = obj.y(indMin:indMax);
			
			method = 'Kasa';
			switch method
				case 'Kasa'
					[xc,yc,R,e] = fitCircle_Kasa(xsub, ysub);
				otherwise
					% 
			end%switch
			
			% Create WAYPOINTS object
			phi = linspace(0, 2*pi, N);
			circle = Waypoints(...
				R*cos(phi) + xc, ...
				R*sin(phi) + yc, ...
				2*R*phi, ...
				1/R*ones(1,N), ...
				phi + pi/2, ...
				1, 1);
			
			% Plot if requested
			if doPlot
				plot(obj, 'b-', 'MarkerSize',10, 'DisplayName','Waypoints');
				hold on
				plot(xsub, ysub, 'r.', ...
					'DisplayName','Waypoints under test');
				% Plot fitted circle and every 10th tangent
% 				plottangent(circle, linspace(1, N, (N-1)/10+1), ...
% 					'Color','k', 'DisplayName','Circle');
				plot(circle, 'Color','k', 'DisplayName','Circle');
				hold off
			end%if
			
		end%fcn
		
		function obj = getSamples(obj, idx)
		%GETSAMPLES		Select subset of waypoints.
		%	OBJ = GETSAMPLES(OBJ,IDX) selects a subset of WAYPOINTS object
		%	OBJ specified by indices IDX.
		%	
		%	Indices IDX must be monotonically increasing!
		%
		%	See also WAYPOINTS/SETSTARTINDEX, WAYPOINTS/NBR.
			
			%%% handle input arguments
			narginchk(2, 2);
			
			if isempty(idx)
				return
			end%if
			if any(~isfinite(idx))
				error('WAYPOINTS:getSamples', 'Indexes must be finite!');
			end%if
			
			if max(idx(:)) > numwp(obj)
				error('WAYPOINTS:getSamples',...
					'Max. index (%d) exceeds number of waypoints(%d).',...
					max(idx(:)), numwp(obj));
			end%if
			
			if any(diff(idx(:)) < 0)
				error('WAYPOINTS:getSamples',...
					'Indices IDX must be monotonically increasing!');
			end%if
			
			
			%%% Get subset of waypoints
			obj = Waypoints(...
				obj.x(idx), ...
				obj.y(idx), ...
				obj.s(idx) - obj.s(idx(1)), ...
				obj.Curv(idx), ...
				obj.Head(idx), ...
				obj.Type(idx), ...
				obj.Nbr(idx) - obj.Nbr(idx(1)) + 1);
			
		end%fcn

		function N = numwp(obj)
		%NUMWP	Number of waypoints.
		%	N = NUMWP(OBJ) returns the number N of waypoints of object OBJ. 
			N = numel(obj.x);
		end%fcn
		
		function d = perpendicularDistance(obj, P1, P2, doPlot, indPlot)
		%PERPENDICULARDISTANCE	Perpendicular distance to line.
		%	D = PERPENDICULARDISTANCE(OBJ,P1,P2) calculate the
		%	perpendicular distance D for all waypoints of OBJ to the line
		%	passing through P1 and P2.
		% 
		%	D = PERPENDICULARDISTANCE(...,DOPLOT) also shows a plot of
		%	results if DOPLOT evaluates to TRUE.
		%
		%	NOTE: The distance for waypoints of OBJ to the left/right of
		%	the line from P1 to P2 is positive/negative.
		%	
		%	See also
		%	https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Line_defined_by_two_points
			
			% Handle input arguments
			narginchk(1, 5)
			
			if nargin < 2
				P1 = obj.InitPoint;
			end%if
			if nargin < 3
				P2 = obj.TermPoint;
			end%if
			
			if (numel(P1) ~= 2) || (numel(P2) ~= 2)
				error('WAYPOINTS:perpendicularDistance',...
					'Line must be defined using two 2-D points.');
			elseif all(P1(:) == P2(:))
				error('WAYPOINTS:perpendicularDistance',...
					'Points are equal!');
			end%if
			
			if nargin < 4
				doPlot = false;
			end%if
			
			
			% Calculate the perpendicular distance for all waypoints of
			% OBJ. See:
			% https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Line_defined_by_two_points
			% Scale by -1 so that points to the left/right of the line from
			% P1 to P2 have a positive/negative distance value.
			d = -((P2(2)-P1(2))*obj.x - (P2(1)-P1(1))*obj.y + P2(1)*P1(2) - P2(2)*P1(1)) ...
				/ sqrt((P2(2)-P1(2))^2 + (P2(1)-P1(1))^2);
			
			
			%%% Visualization
			if doPlot
				if nargin < 5
					indPlot = 1:numwp(obj);
				end%if
				if max(indPlot) > numwp(obj)
					warning('Plot indices limited to the number %d of waypoints!', numwp(obj));
					indPlot = indPlot(indPlot <= numwp(obj));
				end%if
				
				% plot the waypoints
				plot(obj, 'b.', 'MarkerSize',10);
				hold all
				
				% plot the line defined by points P1/P2
				plot([P1(1) P2(1)], [P1(2) P2(2)], 'k', ...
					'LineWidth',1, 'Marker','o');
				
				% plot the perpendicular lines from waypoints to line P12
				% defined by points P1/P2; first calculate its slope and
				% offset
				k12 = (P2(2) - P1(2))/(P2(1) - P1(1));
				d12 = P1(2) - k12*P1(1);
				
				if isfinite(k12) && (k12 ~= 0)
					dstar = obj.y(indPlot) + 1/k12*obj.x(indPlot);
					x_ = (dstar-d12)/(k12+1/k12);
					y_ = k12*x_ + d12;
				elseif ~isfinite(k12)
					% x-coordinates of P1 and P2 are the same
					x_ = P1(1)*ones(size(obj.x(indPlot)));
					y_ = 0*x_ + obj.y(indPlot);
				else
					% slope of line P12 is 0
					x_ = obj.x(indPlot);
					y_ = k12*x_ + d12;
				end
				plot([obj.x(indPlot); x_], [obj.y(indPlot); y_], 'r.-');
				
				hold off
			end%if
			
		end%fcn
		
		function obj = plus(obj1, obj2)
		%+ Plus.
		%	OBJ = OBJ1 + OBJ2 appends waypoints OBJ2 to waypoints OBJ1 in a
		%	way that the terminal point of OBJ2 and initial point of OBJ1
		%	match, keeping only one of these two waypoints.
		%	
		%	Therefore, OBJ consists of one waypoint less than the total
		%	number of waypoints of OBJ1 and OBJ2.
		%	
		%	Note that here plus (+) is a non-commutative operation!
		
			if numwp(obj2) < 2
				obj = obj1;
				return
			end%if
			
			obj = Waypoints(...
				[obj1.x;	obj1.x(end) + obj2.x(2:end) - obj2.x(1)],...
				[obj1.y;	obj1.y(end) + obj2.y(2:end) - obj2.y(1)],...
				[obj1.s;	obj1.s(end) + obj2.s(2:end)],...
				[obj1.Curv;	obj2.Curv(2:end)],...
				[obj1.Head;	obj2.Head(2:end)],...
				[obj1.Type; obj2.Type(2:end)],...
				[obj1.Nbr;	obj2.Nbr(2:end) - obj2.Nbr(1) + 1 + obj1.Nbr(end)]);
			
		end%fcn
		
		function obj = resample(obj, ds, method)
        %RESAMPLE   Resample waypoints.
        %   OBJ = RESAMPLE(OBJ,DS) resamples waypoints OBJ with desired
        %   maximum distance DS between consecutive waypoints by
        %   interpolating x/y data.
		%	
		%	OBJ = OBJ(...,METHOD) specifies the interpolation method.
		%	Default is 'linear'.
		% 
		%	Notes:
		%	 - Length may increase/decrease due to finer/coarser resampling
		%	 for non-straight waypoints.
		% 
		%	See also INTERP1.
			
			% Check input data
			if numwp(obj) < 2
				error('Resampling requires at least two waypoints!');
			end%if
			
			% Handle input arguments
			narginchk(2,3);
			if nargin < 3
				method = 'linear';
			end%if
			
			% Validate input arguments
			validateattributes(ds, {'numeric'}, {'scalar', 'finite'});
			validateattributes(method, {'char'}, {'vector'});
			
			% Return with unmodified WAYPOINTS object
			if ds <= 0
				return
			end%if 
						
			% Number of points required to match DS
			nbrPoints_set = ceil(obj.s(end)/ds) + 1;
			
			% Consider waypoints as parameterized with curve length as the
			% independent parameter
			t = obj.s;
			for i = 1:100 % run at max 100 iterations
				% get interpolated x/y values at query points TQ
				tq = linspace(t(1), t(end), nbrPoints_set)';
				XY = interp1(t, [obj.x, obj.y], tq, method);
% 				S = [0, cumsum(sqrt(diff(XY(1,:)).^2 + diff(XY(2,:)).^2))];
				S = sFrom_x_y(XY(:,1), XY(:,2));
				
				% due to finer/coarser resampling, the length
				% increases/decreases and therefore DELTAACT <= DS might
				% not be fulfilled
				deltaAct = S(end)/(nbrPoints_set - 1);
				
				if deltaAct > ds
					nbrPoints_set = nbrPoints_set + 1;
				else
					break;
				end%if
			end%for
			
			if deltaAct > ds
				fprintf('Maximum number of %u iterations not sufficient to match DS.\n', ...
					uint16(i));
			end%if
			
			% Try interpolating curvature and orientation. Depending on the
			% method, this might fail if the function values contain NaN's
			% TODO: NaNs should are not allowed at constructor
			K = interp1(t, obj.Curv, tq, method);
			PHI = interp1(t, obj.Head, tq, method);
			
			% Resample curvature Type
			TYPE = resample_on_s(obj, 'Type', S);
			
			% Resample segment number
			NBR	= resample_on_s(obj, 'Nbr', S);
			
			% Create the resampled WAYPOINTS object
			obj = Waypoints(XY(:,1), XY(:,2), S, K, PHI, TYPE, NBR);
			
		end%fcn
		
		function obj = reverse(obj)
		%REVERSE	Reverse order of waypoints.
		%	OBJ = REVERSE(OBJ) reverses the direction of waypoints OBJ so
		%	[x(end),y(end)] becomes [x(1),y(1)] and so on.
		%	
		%	NOTE: This method supports array inputs OBJ!
			
			% Handle input arguments
			narginchk(1,1);
			
			% Reverse waypoints
			for i = 1:numel(obj)
				% Cast NBR property to sign preserving data Type!
				nbrS = double(obj(i).Nbr);
				
				obj(i) = Waypoints(...
					+flip(obj(i).x),...
					+flip(obj(i).y),...
					-flip(obj(i).s) + obj(i).s(end),... % in case of non-equally distributed x/y
					-flip(obj(i).Curv),...
					+flip(obj(i).Head) + pi,... 
					+flip(obj(i).Type),...
					-flip(nbrS) + nbrS(end) + 1); 
			end%for
			
		end%fcn
		
		function obj = rotate(obj, phi)
		%ROTATE		Rotate WAYPOINTS object.
		%	OBJ = ROTATE(OBJ,PHI) rotates OBJ by an angle PHI in radians
		%	around the origin.
		% 
		%	OBJ = ROTATE(OBJ) applies the default value -OBJ(1).Head(1) for
		%	PHI.
		%	
		%	If OBJ is an array of WAYPOINTS objects, PHI is applied to all
		%	WAYPOINTS object.
			
			% Handle input arguments
			narginchk(1, 2);
			
			if nargin < 2
				phi = -obj(1).Head(1);
			end%if
			
			if numel(phi) ~= 1 || ~isnumeric(phi)
				error(['Method ROTATE requires a numeric input',...
					' argument with one element.']);
			end%if
			
			
			% Rotate
			% rotation matrix in R^2.
			rotMat = [cos(phi) -sin(phi); 
					  sin(phi) +cos(phi)]';
			
			for i = 1:numel(obj)
				% much faster than rotMat*[obj.x,obj.y]'
				xy_new = [obj(i).x, obj(i).y]*rotMat;
				
				obj(i) = Waypoints(...
					xy_new(:, 1), ...
					xy_new(:, 2), ...
					obj(i).s, ...
					obj(i).Curv, ...
					obj(i).Head + phi, ...
					obj(i).Type, ...
					obj(i).Nbr);
			end%for
			
		end%fcn
		
		function obj = setStartIndex(obj,idx)
		%SETSTARTINDEX	Set waypoint as initial waypoint.
		%	OBJ = SETSTARTINDEX(OBJ,INDX) reorders waypoints so that index
		%	INDX becomes the inital waypoint.
		%	
		%	Property NBR is manipulated according to it's requirements!
		%	
		%	WARNING: this method only makes sense for waypoints
		%	representing closed paths circuits!
		% 
		%	See also WAYPOINTS/NBR.
			
			%%% Handle input arguments
			if idx == 1
				return;
			end%
			
			% Distance from last to first element
			delta = sqrt((obj.x(1)-obj.x(end))^2 + (obj.y(1)-obj.y(end))^2);
			
			% The length should be stricly monotonically increasing
			s_new = [...
				obj.s(idx:end) - obj.s(idx);...
				obj.s(1:idx-1) + obj.s(end) - obj.s(idx) + delta];
			
			%%% Reorder waypoints
			% Make sure that property NBR is monotonically increasing and
			% starts at 1
			obj = Waypoints(...
				[obj.x(idx:end);	obj.x(1:idx-1)],...
				[obj.y(idx:end);	obj.y(1:idx-1)],...
				s_new,...
				[obj.Curv(idx:end);	obj.Curv(1:idx-1)],...
				[obj.Head(idx:end);	obj.Head(1:idx-1)],...
				[obj.Type(idx:end); obj.Type(1:idx-1)],...
				[obj.Nbr(idx:end);	obj.Nbr(1:idx-1) + obj.Nbr(end)] - obj.Nbr(idx) + 1);
			
		end%fcn
		
		function obj = shiftBy(obj, P)
		%SHIFTBY	Shift waypoints by offset.
		%	OBJ = SHIFTBY(OBJ,P) shifts waypoints so that the initial point
		%	is [OBJ.x(1)+P(1) OBJ.y(1)+P(2)].
		% 
		%	NOTE: If OBJ is an array of WAYPOINTS, P is applied to all
		%	array elements!
			
			% Handle input arguments
			narginchk(2, 2);
			
			if numel(P) ~= 2 || ~isnumeric(P)
				error(['Method SHIFTBY requires a numeric input',...
					' argument with two elements.']);
			end%if
			
			% Shift waypoints
			for i = 1:numel(obj)
				obj(i) = Waypoints(...
					obj(i).x + P(1),...
					obj(i).y + P(2),...
					obj(i).s,...
					obj(i).Curv,...
					obj(i).Head,...
					obj(i).Type,...
					obj(i).Nbr);
			end%for
			
		end%fcn
		
		function obj = shiftTo(obj, P)
		%SHIFTTO	Shift waypoints to point.
		%	OBJ = SHIFTTO(OBJ,P) shifts waypoints to new inital point P.
		%	
		%	OBJ = SHIFTTO(OBJ) applies the default value [0 0] for P.
		%	
		%	NOTE: If OBJ is an array of WAYPOINTS, P is applied to the
		%	first array element, the others are shifted to maintain their
		%	positioning with respect to the first WAYPOINTS object!
			
			
			% Handle input arguments
			narginchk(1, 2);
			
			if nargin < 2
				P = [0 0];
			end%if
			
			if numel(P) ~= 2 || ~isnumeric(P)
				error(['Method SHIFTTO requires a numeric input',...
					' argument with two elements.']);
			end%if
			
			% Shift waypoints and maintain positional relation if OBJ is an
			% array
			P0 = obj(1).InitPoint;
			for i = 1:numel(obj)
				Pi = P + obj(i).InitPoint - P0;
				obj(i) = Waypoints(...
					obj(i).x - obj(i).x(1) + Pi(1),...
					obj(i).y - obj(i).y(1) + Pi(2),...
					obj(i).s,...
					obj(i).Curv,...
					obj(i).Head,...
					obj(i).Type,...
					obj(i).Nbr);
			end%for
			
		end%fcn
		
		function obj = unwrap(obj)
		%UNWRAP		Unwrap heading angle property HEAD.
		%	OBJ = UNWRAP(OBJ) applies UNWRAP() to property HEAD.
		%	
		%	NOTE: This method supports array inputs OBJ!
		% 
		%	See also UNWRAP.
			
			for i = 1:numel(obj)
				obj(i).Head = unwrap( obj(i).Head );
			end%for
			
		end%fcn
		
		function write2file(obj, fn, format, varargin)
		%WRITE2FILE		Write waypoints to file.
		%	WRITE2FILE(OBJ,FN) writes waypoints OBJ to file with filename
		%	FN (specify extension!).
		%	
		%	WRITE2FILE(OBJ,FN,FORMAT,VARARGIN) lets you specify the output
		%	format. Supported values for FORMAT are:
		%	  - 'raw' (default value)
		%		 One column per property of OBJ. Property names are used as
		%		 labels in the first row.
		%		 
		%	  - 'CarMaker'
		%		Depending on the file name extension: 
		%		 *.csv: ASCII file (cartesian coordinates), supports at
		%		 least CarMaker v4.0 - v5.0.1. Use VARARGIN to set the left
		%		 and right lane widths/margin widths: [lw_left lw_right],
		%		 [mw_left mw_right]. If unspecified, values are set to
		%		 zero.
		%		 
		%		 *.kml: KML file (WGS84-coordinates), supports at least
		%		 CarMaker v4.0 - v5.0.1. Use VARARGIN to specify a scalar
		%		 or vector of UTM zone integers. 
		%
			
		
			% Set the default FORMAT
			if nargin < 3 || isempty(format)
				format = 'raw';
			end%if
			
			% Open file
			fid = openFileWithExt(fn);
			
			switch format
				case 'raw'
					fprintf(fid,...
						'%12s %12s %12s %12s %12s %4s %3s \n',...
						'x [m]','y [m]','s [m]','k [1/m]','phi [rad]','type','nbr');
					fprintf(fid,...
						'%+12.3f %+12.3f %12.3f %12.3f %12.3f %4d %3u \n',...
						[obj.x; obj.y; obj.s; obj.Curv; obj.Head;...
						double(obj.Type); double(obj.Nbr)]);
					
				case 'CarMaker'
					write2file_CarMaker(fn, obj, varargin{:});
					
				case 'rompac'
					fprintf(fid,...
						'%s %10s %12s %12s %12s %12s \n',...
						'#', 'x_[m]', 'y_[m]', 'phi_[rad]','SD_right_[m]','SD_left_[m]');
					fprintf(fid,...
						'%+12.3f %+12.3f %12.3f %12.3f %12.3f \n',...
						[obj.x; obj.y; obj.Head; varargin{1}']);
					
				otherwise
					error('WAYPOINTS:write2File','Unknown FORMAT specifier');
					
			end%switch
			
			% Close file
			fclose(fid);
			
		end%fcn
		
		function [out,latOff_LAD,angDev_LAD,curvat_LAD,isValid_LAD] = ...
				laneTracking(obj,xyCG_global,yawAngle_global,LAD,mode)
		%LANETRACKING	Calculate lane tracking pose.
		%	
		%	[~,LATOFF_LAD,ANGDEV_LAD,CURVAT_LAD,ISVALID_LAD] =
		%	LANETRACKING(OBJ,XYCG_GLOBAL,YAWANGLE_GLOBAL,LAD) calculates
		%	the lateral offset LATOFF_LAD, the angular deviation ANGDEV_LAD
		%	and the curvature CURVAT_LAD at the look-ahead distances LAD in
		%	front of the vehicles center of gravity XYCG_GLOBAL along its
		%	longitudinal axis oriented with the angle YAWANGLE_GLOBAL and
		%	with respect to waypoints OBJ. The logical vector ISVALID_LAD
		%	indicates if the corresponding entry of all other output
		%	arguments is valid (true) or not (false).
		%	
		%	OUT = LANETRACKING(...) is the syntax to be used from simulink,
		%	where OUT = [LAD;ISVALID_LAD;LATOFF_LAD;ANGDEV_LAD;CURVAT_LAD].
		%	
		%	The number of columns of all output arguments matches the
		%	length of input argument LAD, moreover the i-th column of any
		%	output argument corresponds to the i-th element of look-ahead
		%	distance LAD.
		
		% Subject: lka
			
			% Methode: Koordinatentransformation
			if nargin <5
				mode = 0;
			end%if
			
			% ensure row/column format
			LAD = LAD(:)';
			xyCG_global = xyCG_global(:);
			
			%%% shift origin to vehicles CG/rotate so vehicle is oriented
			%%% along global x-axis
			% CG (transformed)
			xyCG_T = [0;0]; % xyCG_global - xyCG_global
			
			% Sollbahn (transformed)
			obj_T = shiftTo(obj, [obj.x(1);obj.y(1)]-xyCG_global);
			obj_T = rotate(obj_T, -yawAngle_global);
			
			% Koordinaten des Punkts bei look-ahead distance (transformed)
			xyLAD_T = [xyCG_T(1) + LAD; xyCG_T(2)+zeros(size(LAD))];
			
			
% 			%%% shift origin to point at LAD
% 			xyCG_T	= xyCG_T - [lad;0];
% 			obj_T	= shiftTo(obj_T,[obj_T.x(1);obj_T.y(1)] - [lad;0]);
% 			xyLAD_T	= xyLAD_T - [lad;0];
			
			
			%%% Indexes of waypoints relevant for tracking error
			%%% calculation
			% maximaler Abstand zwischen x-Werten der Sollbahn
			deltaX_max = max(abs(diff(obj_T.x)));
			
			% Depending on the shape of the set of waypoints (e.g. closed
			% path), there might be multiple waypoints whose x-coordinates
			% are within the range of LAD +- DELTAX_MAX/2.
			% 
			% Get logical indices where OBJ_T.X is in the range of LADs
			% x-coordinate
			%  rows: LAD 
			%  columns: waypoint x-coordinates
			logIndx1 = bsxfun(@le, xyLAD_T(1,:)' - deltaX_max/2, obj_T.x');
			logIndx2 = bsxfun(@ge, xyLAD_T(1,:)' + deltaX_max/2, obj_T.x');
			logIndx = logIndx1 & logIndx2;
			
			% error if no element of LOGINDX is true
			if ~any(any(logIndx))
				plotLaneTracking(obj,xyCG_global,yawAngle_global,LAD,[],obj_T,xyCG_T);
				error('WAYPOINTS:laneTracking',...
					['Keine Elemente der Sollbahn im Bereich der',... 
					' aktuellen Fahrzeugposition gefunden'])
			end%if
			
			% Numerical indices according to LOGINDX: NUMINDROW points to
			% the according LAD, NUMINDCOL points to the according
			% x-coordinate (see LOGINDX calculation above)
			[numIndRow,numIndCol] = find(logIndx);
			
			if mode
				%%% reduce number of potential elements
				% since there might be multiple elements within LAD +-
				% DELTAX_MAX/2, get the closest in x-direction per LAD AND
				% per contigous indices
				
				% split NUMINDCOL into segments of contigous indices and
				% LAD values
				splitByLAD = find(diff(numIndRow)>0);
				splitByInd = find(diff(numIndCol)>1);	
				splitInd = [0;unique([splitByLAD;splitByInd]);length(numIndRow)];
				
				% get index of element closest to LAD, do it for each
				% index-segment defined by SPLITIND
				numIndRow_red = zeros(length(splitInd)-1,1);
				numIndCol_red = zeros(length(splitInd)-1,1);
				for i = 1:length(splitInd)-1
					numIndRow_i = numIndRow(splitInd(i)+1:splitInd(i+1));
					numIndCol_i = numIndCol(splitInd(i)+1:splitInd(i+1));
					
					% since we grouped by contigous indices and LAD,
					% NUMINDWOR_I must not contain different values
					if ~all(numIndRow_i(1) == numIndRow_i)
						error('This should not have happened!');
					end%if
					
					%
					if length(numIndRow_i) > 2
						[~,winInd] = min(abs(obj_T.x(numIndCol_i) - xyLAD_T(1,numIndRow_i(1))));
					else
						winInd = 1;
					end%if
					
					numIndRow_red(i) = numIndRow_i(winInd);
					numIndCol_red(i) = numIndCol_i(winInd);
				end%for
				numIndRow = numIndRow_red;
				numIndCol = numIndCol_red;
			end%if
			
% 			plotLaneTracking(obj,xyCG_global,yawAngle_global,LAD,numIndCol,obj_T,xyCG_T);
			
			
			%%% get lateral offset for potential elements
			% interpolation index-range: m>1 mainly increase the
			% calculation time, lateral offset and angular deviation are
			% rarely affected.
			m = 1;
			[indl,indu] = interpIndexRange(numIndCol, [1,numwp(obj)], m);
			
			% preallocation of for-loop variable
			lanePose_LAD_candidates = zeros(3,length(numIndCol));
			try
				% Inter- bzw. Extrapoliere y-Werte der transf. Solltrajektorie
				for i = 1:length(numIndCol)
					% same result like interp1(..,'spline') but faster
					lanePose_LAD_candidates(:,i) = spline(...
						obj_T.x(indl(i):indu(i)),...
						[obj_T.y(indl(i):indu(i)),...
						 obj_T.Head(indl(i):indu(i)),...
						 obj_T.Curv(indl(i):indu(i))],...
						xyLAD_T(1, numIndRow(i)));
				end%for
			catch exception
				plotLaneTracking(obj, xyCG_global, yawAngle_global, ...
					LAD, numIndCol, obj_T, xyCG_T);
				error(exception.message);
			end%try
			latOff_LAD_candidates = lanePose_LAD_candidates(1,:);
			angDev_LAD_candidates = lanePose_LAD_candidates(2,:);
			curvat_LAD_candidates = lanePose_LAD_candidates(3,:);
			
			
			%%% select one of multiple lateral offsets per LAD-value
			% get most possible (by means of smallest) value per LAD-value
			latOff_LAD	= zeros(size(LAD));
			angDev_LAD	= zeros(size(LAD));
			curvat_LAD	= zeros(size(LAD));
			isValid_LAD = false(size(LAD));
			for i = 1:length(latOff_LAD_candidates)
				latOff_underTest = latOff_LAD_candidates(i);
				angDev_underTest = angDev_LAD_candidates(i);
				curvat_underTest = curvat_LAD_candidates(i);
				LAD_group = numIndRow(i);
				
				if abs(latOff_underTest) < abs(latOff_LAD(LAD_group)) || ...
						~isValid_LAD(LAD_group)
					latOff_LAD(LAD_group)	= latOff_underTest;
					angDev_LAD(LAD_group)	= angDev_underTest;
					curvat_LAD(LAD_group)	= curvat_underTest;
					isValid_LAD(LAD_group)	= true;
				end%if
				
			end%for
			
			% set multiple occurences of same LAD values to invalid
			[LAD_unique,ia] = unique(LAD,'stable');
			if numel(LAD_unique) ~= numel(LAD)
				% duplicate indices
				duplicate_ind = setdiff(1:numel(LAD),ia);
				
				% set to false
				isValid_LAD(duplicate_ind) = false;
			end%if
			
			% collect all output arguments (to be used in simulink)
			out = [LAD;isValid_LAD;latOff_LAD;angDev_LAD;curvat_LAD];
			
		end%fcn
		
	end%methods
	
	
	methods (Access = private)
		
		[h,axh] = plot_raw(axh, obj, varargin)
		[h,axh] = quiver_raw(axh, obj, varargin)
		PROP = resample_on_s(obj, field, S)
		cellStr = getTitleCellString(obj)
		write2file_CarMaker(fn, obj, varargin)
		
	end%methods
	
	
	%%% GET-Methods
	methods
		
		function pos = get.InitPoint(obj)
			pos = [obj.x(1); obj.y(1)];
		end%fcn
		
		function pos = get.TermPoint(obj)
			pos = [obj.x(end); obj.y(end)];
		end%fcn
		
	end%methods
	
	
	%%% SET-Methods
	methods
		
		function obj = set.Name(obj, value)
			if ischar(value)
				obj.Name = value;
			else
				error('Name argument must be of class char!')
			end
		end%fcn
		
		function obj = set.s(obj, value)
			% The length can not decrease from one to following point
			if any(diff(value) < 0)
				error('WAYPOINTS:DecreasingS', ...
					'Property S must be monotonically increasing!')
			end%if
			
			obj.s = value;
		end%fcn
		
		function obj = set.Type(obj, value)
			
			% Check value range
			if any(value < -1) || any(value > 4)
				error(['Property TYPE out of range! ',...
					'Type help WAYPOINTS/TYPE for valid values.']);
			end%if
			
			% Ensure class
			if ~isa(value, 'int8')
				obj.Type = int8(value);
			else
				obj.Type = value;
			end%if
			
		end%fcn
		
		function obj = set.Nbr(obj, value)
			
			% NBR must be monotonically increasing
			assert(~any(diff(value) < 0), ...
				'Input NBR must be monotonically increasing!');
			
			% Make sure NBR starts with 1
			if value(1) ~= 1
				fprintf('Warning: Modifying property NBR to start at 1!');
				value = value - value(1) + 1;
			end%if
			
			% Ensure class
			if ~isa(value, 'uint16')
				obj.Nbr = uint16(value);
			else
				obj.Nbr = value;
			end%if
			
		end%fcn
		
	end%SET-Methods
	
	
	%%% Static-Methods
	methods (Static)
		
		function obj = ll2Waypoints(lat, lon, name)%#codegen
		%LL2WAYPOINTS	Convert LAT/LON coordinates to WAYPOINTS.
		%	OBJ = LL2WAYPOINTS(LAT,LON)
			
			% LAT/LON are required, accept 3 optional name-value pair
			% arguments
			narginchk(2, 3);
			
			% Convert from lat/lon to UTM
			[x,y] = ll2utm(lat(:), lon(:));
			s = sFrom_x_y(x, y);
			
			head = cx2Heading(gradient(x), gradient(y));
			curv = curvFrom_head_s(head, s);
			
			% Create WAYPOINTS object
			obj = Waypoints(x, y, s, curv, head);
			if nargin > 2
				obj.Name = name;
			end
			
		end%fcn
		
		function obj = pp2Waypoints(pp, t, name, mode)
		% PP2WAYPOINTS	Convert piecewise polynomial structure to WAYPOINTS. 
		%   OBJ = PP2WAYPOINTS(PP,T) creates WAYPOINTS instace OBJ from
		%   piecewise polynomial structures PP sampled at T.
		% 
		%	OBJ = PP2WAYPOINTS(PP,T,NAME,MODE) setting MODE='integral,
		%	length S is calculated via numerical integration. Default value
		%	is 'cumsum'.
		%	
		%	See also MKPP, UNMKPP, SPLINE.
			
			narginchk(2, 4);
			
			if nargin < 4
				mode = 'cumsum';
			end%if
			
			% Get dimension via UNMKPP to check for valid piecewise
			% polynomial struct.
			[~,~,~,~,ppDim] = unmkpp(pp);
			% PP needs to return a 2-D array (x and y).
			if ppDim ~= 2
				error('Piecewise polynomial PP must have a dimension of 2!');
			end%if
			
			% Make sure T is a row vector. As a result, PPVAL will return a
			% matrix whose number of columns matches the number of elements
			% of T.
			if ~isvector(t) || ~isrow(t)
				error('Query point must be a row vector!');
			end%if
			
			% Evaluate piecewise polynomial structs
			ppd1 = ppdiff(pp, 1);
			ppd2 = ppdiff(ppd1, 1);
			% PPVAL returns an array of size DIM-by-numel(T)
			xy	 = ppval(pp, t)';
			d1xy = ppval(ppd1, t)';
			d2xy = ppval(ppd2, t)';
			
			switch mode
				case 'cumsum'
					s = sFrom_x_y(xy(:,1), xy(:,2));
					
				case 'integral'
					% Numerical integration of path length
					fun = @(t) sqrt(sum(ppval(ppd1, t).^2, 1));
					s = zeros(numel(t), 1);
					for i = 1:numel(t)-1
						s(i+1) = integral(@(t) fun(t), t(i), t(i+1));
					end%for
					s = cumsum(s);
					
				otherwise
					error('Unknown mode!')
			end%%switch
			
			head = cx2Heading(d1xy(:,1), d1xy(:,2));
			curv = cx2Curvature(d1xy(:,1), d1xy(:,2), d2xy(:,1), d2xy(:,2));
			
			obj = Waypoints(xy(:,1), xy(:,2), s, curv, head);
			if nargin > 2
				obj.Name = name;
			end
			
		end%fcn
		
		function obj = xy2Waypoints(x, y, name)
		%XY2WAYPOINTS	Convert x/y coordinates to WAYPOINTS. 
		%   OBJ = XY2WAYPOINTS(X,Y) creates WAYPOINTS instace OBJ from
		%   cartesian coordinates X and Y, which must be vectors with the
		%   same number of elements!
		%
		%	Property S is calculated using cumulative sum.
		%	Properties PHI and CURV are calculated as difference quotients.
		%	Property TYPE is set to -1.
		%	Property NBR is set to 1.
		
		
			narginchk(2, 3);
		
			%%% Handle input arguments
			if ~isvector(x) || ~isvector(y)
				error('Input arguments X/Y must be vectors!');
			end%if
			
			if numel(x) ~= numel(y)
				error('Input arguments X/Y must have the same number of elements!');
			else
				x = x(:);
				y = y(:);
			end%if
			
			
			% Check for NaN's: GRADIENT will fail to return a numeric value
			% for terminal points in case of NANs
			isNan = or(isnan(x), isnan(y));
			
			%%% Convert to WAYPOINTS
			s = sFrom_x_y(x(~isNan), y(~isNan));
			head = cx2Heading(gradient(x(~isNan)), gradient(y(~isNan)));
			curv = curvFrom_head_s(head, s);
			
			obj = Waypoints(x(~isNan), y(~isNan), s, curv, head, -1, 1);
			if nargin > 2
				obj.Name = name;
			end
			
		end%fcn
		
		function obj = circle(r, phi01, N)
		%CIRCLE		Create circle.
		%	OBJ = WAYPOINTS.CIRCLE(R) creates WAYPOINTS instance OBJ
		%	describing a circle of radius R.
		%
		%	OBJ = WAYPOINTS.CIRCLE(R, PHI01) sets the inital and final
		%	angle to PHI01(1) and PHI01(2) respectively. Defalut value is
		%	[0; 2*pi];
		%
		%	OBJ = WAYPOINTS.CIRCLE(R, PHI01, N) creates the circle using N
		%	samples. Default value is N = 100.
			
			if nargin < 3
				N = 100;
			end
			if nargin < 2
				phi01 = [0; 2*pi];
			end
			t = linspace(phi01(1), phi01(2), N);
			obj = Waypoints(r*cos(t), r*sin(t), r*(t-t(1)), repmat(1/r,N,1), t+pi/2);
			obj.Name = ['Circle R=',num2str(r)]; 
		end%fcn
		
	end%methods
	
end%class
