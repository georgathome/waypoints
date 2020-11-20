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
%	 CURV	- Curvature [1/m].
%	 HEAD	- Heading angle [rad].
%	 TYPE	- Integer indicating the waypoint type.
%	 NBR	- Segment number of joined WAYPOINTS.
%	 INITPOINT - Initial point [x(1) y(1)].
%	 TERMPOINT - Terminal Point [x(end) y(end)].
%	 NBREL	- Number of waypoints in the WAYPOINTS object.
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
%	   douglasPeuker - Ramer�Douglas�Peucker point reduction algorithm. 
%	   fitStraight	 - Fit straight line to waypoints.
%	   fitCircle	 - Fit circle to waypoints.
%	 
%	 - VISUALIZATION
%	   curvplot	 - Plot waypoints, heading and curvature.
%	   plot		 - Plot waypoints.
%	   plotdiff	 - Plot waypoints with specific appearance.
%	   plotdiff_ - Plot waypoints properties with specific appearance.
%	   plottangent - Plot waypoints and tangents.
%	   quiver	 - ToDo
%	 
%	 - EXPORT
%	   write2file - Write WAYPOINTS properties to file.
% 
%	 - MISC
%	   laneTracking - Get the lane tracking pose.
%	 
%	
%	See also LKASEGMENT, LKASEGMENTSTRAIGHT, LKASEGMENTCIRCLE,
%	LKASEGMENTCLOTHOID.
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



	properties (Constant, Hidden)
		% Curvature types
		curvTypes = {'waypoints','straight','circle','clothoid','sine'};
	end%properties
	
	properties
		% NAME - The name of the WAYPOINTS object (char or string).
		Name = '';
	end%properties
	
	properties (SetAccess = private)
		% X - 1-by-n double of X-coordinates [m].
		x = [0 1];
		
		% Y - 1-by-n double of Y-coordinates [m].
		y = [0 0];
		
		% S - 1-by-n double of lengths [m].
		s = [0 1];
		
		% CURV - 1-by-n double of curvatures [1/m].
		curv = [0 0];
		
		% HEAD - 1-by-n double of heading angles [rad].
		%	 Measured from the positive x-axis to the waypoint's tangent.
		head = zeros(1, 2);
		
		% TYPE - 1-by-n int8 of curvature types [-]:
		%	-1 .. unknown/undefined
		%	 0 .. straight
		%	 1 .. circular
		%	 2 .. clothoid
		type = zeros(1, 2, 'int8');
		% See also WAYPOINTS/CURVTYPES
		
		% NBR - 1-by-n uint16 of segment number [-].
		%	When appending multiple WAYPOINTS objects to one single
		%	WAYPOINTS object, this property allows to distinguish them in
		%	the single object.
		%
		%	Needs to be monotonically increasing and start with 1!
		nbr = ones(1, 2, 'uint16');
	end%properties
	
	properties (Dependent)
		% INITPOINT - The initial waypoint.
		initPoint
		
		% TERMPOINT - The terminal (end) waypoint.
		termPoint
		
		% NBREL - The number of waypoints the WAYPOINTS object consists of.
		nbrel
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
			sz = size(x);
			if (sz(1) ~= 1) || (sz(2) < 1)
				error('WAYPOINTS:Waypoints:nonemptyColumnVectors',...
					'Inputs have to be non-empty row vectors!')
			end%if
			if ~isequal(sz, size(y), size(s), size(head), size(curv))
				error('WAYPOINTS:Waypoints:unequalInputArgumentSizes',...
					'Inputs must have the same size.')
			end%if
			
			% Before we can check the size of optional input arguments TYPE
			% and NBR, we must ensure they are defined
			if nargin < 7
				nbr = ones(sz, 'uint16');
			end%if
			if nargin < 6
				type = -ones(sz, 'int8');
			end%if			
			
			if ~isscalar(type) && ~isequal(sz, size(type))
				error('Property TYPE must be scalar or of size 1-by-%u', sz(2))
			end%if
			if ~isscalar(nbr) && ~isequal(sz, size(nbr))
				error('Property NBR must be scalar or of size 1-by-%u', sz(2))
			end%if
			
			
			%%% (2) Check for NaN's in waypoints
			isNan = or(isnan(x), isnan(y));
			if any(isNan)
				fprintf(['Warning: NaN''s were found in waypoints; ', ...
					'those samples were removed!\n']);
				
				% Update expected size of properties X, ..., CURV
				sz = [1, nnz(~isNan)];
			end%if
			
			
			%%% (3) Set class properties
			obj.x		= x(~isNan);
			obj.y		= y(~isNan);
			obj.s		= s(~isNan);
			obj.head	= head(~isNan);
			obj.curv	= curv(~isNan);
			if isscalar(type)
				% Explicit array size syntax for compatibility in MATLAB
				% Function blocks within Simulink.
				obj.type = type(1) * ones(1, sz(2), class(type));
			else
				obj.type = type(~isNan);
			end%if
			if isscalar(nbr)
				% Same as above!
				obj.nbr = nbr(1) * ones(1, sz(2), class(nbr));
			else
				obj.nbr = nbr(~isNan);
			end%if
			
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
				ds = sqrt( sum((obj2append.initPoint - obj.termPoint).^2) );

				obj = Waypoints(...
					[obj.x,		obj2append.x], ...
					[obj.y,		obj2append.y], ...
					[obj.s,		obj2append.s + obj.s(end) + ds], ...
					[obj.curv,	obj2append.curv], ...
					[obj.head,	obj2append.head], ...
					[obj.type,	obj2append.type], ...
					[obj.nbr,	obj2append.nbr + obj.nbr(end)]);
			end%for
			
		end%fcn
		
		function obj = changeSignOfCurvature(obj)
		% CHANGESIGNOFCURVATURE		Flip waypoints curvature sign.
		%	OBJ = CHANGESIGNOFCURVATURE(OBJ) changes the curvature sign of
		%	waypoints. Properties Y and HEAD are modified accordingly while
		%	maintaining the initial point.
		%	
		%	NOTE: This method supports array inputs OBJ!
		
			for i = 1:numel(obj)
				% Get starting point/angle
				P0 = [obj(i).x(1); obj(i).y(1)];
				phi0 = obj(i).head(1);
			
				% Shift to origin and rotate so initial slope is zero
				obj(i) = shiftTo(obj(i));
				obj(i) = rotate(obj(i), -phi0);
			
				% Invert the curvature
				obj(i) = Waypoints(...
					+obj(i).x,...
					-obj(i).y,...
					+obj(i).s,...
					-obj(i).curv,...
					-obj(i).head,...
					+obj(i).type,...
					+obj(i).nbr);
			
				% Undo the shift/rotate procedure
				obj(i) = rotate(obj(i), phi0);
				obj(i) = shiftTo(obj(i), P0);
			end%for
			
		end%fcn
		
		function [obj,idx] = douglasPeuker(obj, eps)
		% DOUGLASPEUKER		Ramer-Douglas-Peucker point reduction.
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
		%	See also GETSAMPLES.
			
			% Initialize a logical array indicating which waypoints to keep
			N = obj.nbrel;
			keepIdx = true(1, N);
			
			% Recursively set indexes of waypoints that can be discarded to
			% false
			dprec(1, obj.nbrel);
			
			
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
			
			% Handle input arguments
			if (nargin < 2) || isempty(indMinMax)
				indMin = 1;
				indMax = numel(obj.x);
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
			A = [xsub', ones(indMax-indMin+1, 1)];
			b =  ysub';
			
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
				[0, cumsum(sqrt(diff(xsub).^2 + diff(yfit).^2))], ...
				zeros(1, numel(xsub)), ...
				atan2(dy, dx)*ones(1, numel(xsub)), 0, 1);
			
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
			
			% Handle input arguments
			if nargin < 3 || isempty(indMinMax)
				indMin = 1;
				indMax = numel(obj.x);
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
				plot(obj, 'r.', 'MarkerSize',10, 'DisplayName','Waypoints');
				hold on
				plot(xsub, ysub, 'r-', ...
					'DisplayName','Waypoints under test');
				% Plot fitted circle and every 10th tangent
				plottangent(circle, linspace(1, N, (N-1)/10+1), ...
					'Color','k', 'DisplayName','Circle');
				hold off
			end%if			
			
		end%fcn
		
		function obj = getSamples(obj, idx)
		% GETSAMPLES	Select subset of waypoints.
		%	OBJ = GETSAMPLES(OBJ,IDX) selects a subset of WAYPOINTS object
		%	OBJ specified by indices IDX.
		%	
		%	Indices IDX must be monotonically increasing!
		%
		%	See also WAYPOINTS/SETSTARTINDEX, WAYPOINTS/NBR.
			
			%%% handle input arguments
			narginchk(2, 2);
			
			if isempty(idx) || ~isvector(idx) || any(~isfinite(idx))
				error('WAYPOINTS:getSamples',...
					['Input argument INDRANGE must be a non-empty finite vector!',...
					'Type help WAYPOINTS/getSamples.']);
			end%if
			
			if max(idx) > obj.nbrel
				error('WAYPOINTS:getSamples',...
					'Max. index (%d) exceeds number of waypoints(%d).',...
					max(idx), obj.nbrel);
			end%if
			
			if any(diff(idx) < 0)
				error('WAYPOINTS:getSamples',...
					'Indices IDX must be monotonically increasing!');
			end%if
			
			
			%%% Get subset of waypoints
			obj = Waypoints(...
				obj.x(idx), ...
				obj.y(idx), ...
				obj.s(idx) - obj.s(idx(1)), ...
				obj.curv(idx), ...
				obj.head(idx), ...
				obj.type(idx), ...
				obj.nbr(idx) - obj.nbr(idx(1)) + 1);
			
		end%fcn
		
		function dist = lateralDistance(OBJ, obj)
		%UNTITLED Summary of this function goes here
		%   Detailed explanation goes here
		
			N = OBJ.nbrel;
			dist = NaN(N, numel(obj));

			for i = 1:numel(obj)
				for n = 1:N
					[~,latOff_LAD] = laneTracking(obj(i), ...
						[OBJ.x(n); OBJ.y(n)], ... % current position under test
						OBJ.head(n), ... % heading at current position
						0); %LAD
					
					%
					dist(n, i) = latOff_LAD;
				end%for
			end%for
			
		end%fcn
		
		function d = perpendicularDistance(obj, P1, P2, doPlot, indPlot)
		% PERPENDICULARDISTANCE
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
		%	https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Line_defined_by_two_points.
			
			% Handle input arguments
			narginchk(1, 5)
			
			if nargin < 2
				P1 = obj.initPoint;
			end%if
			if nargin < 3
				P2 = obj.termPoint;
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
					indPlot = 1:obj.nbrel;
				end%if
				if max(indPlot) > obj.nbrel
					warning('Plot indices limited to the number %d of waypoints!', obj.nbrel);
					indPlot = indPlot(indPlot <= obj.nbrel);
				end%if
				
				% plot the waypoints
				plot(obj, 'b.', 'MarkerSize',10);
				hold all
				
				% plot the line defined by points P1/P2
				plot([P1(1) P2(1)], [P1(2) P2(2)],'k','LineWidth',1,'Marker','o');
				
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
				plot([obj.x(indPlot); x_], [obj.y(indPlot); y_],'r.-');
				
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
		
			
			if obj2.nbrel < 2
				obj = obj1;
				return
			end%if
			
			obj = Waypoints(...
				[obj1.x,	obj1.x(end) + obj2.x(2:end) - obj2.x(1)],...
				[obj1.y,	obj1.y(end) + obj2.y(2:end) - obj2.y(1)],...
				[obj1.s,	obj1.s(end) + obj2.s(2:end)],...
				[obj1.curv,	obj2.curv(2:end)],...
				[obj1.head,	obj2.head(2:end)],...
				[obj1.type, obj2.type(2:end)],...
				[obj1.nbr,	obj2.nbr(2:end) - obj2.nbr(1) + 1 + obj1.nbr(end)]);
			
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
			if obj.nbrel < 2
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
				tq = linspace(t(1), t(end), nbrPoints_set);
				XY = interp1(t, [obj.x; obj.y]', tq, method)';
				S = [0, cumsum(sqrt(diff(XY(1,:)).^2 + diff(XY(2,:)).^2))];
				
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
			K = interp1(t, obj.curv, tq, method);
			PHI = interp1(t, obj.head, tq, method);
			
			% Resample curvature type
			TYPE = resample_on_s(obj, 'type', S);
			
			% Resample segment number
			NBR	= resample_on_s(obj, 'nbr', S);
			
			% Create the resampled WAYPOINTS object
			obj = Waypoints(XY(1,:), XY(2,:), S, K, PHI, TYPE, NBR);
			
		end%fcn
		
		function obj = reverse(obj)
		% REVERSE	Reverse order of waypoints.
		%	OBJ = REVERSE(OBJ) reverses the direction of waypoints OBJ so
		%	[x(end),y(end)] becomes [x(1),y(1)] and so on.
		%	
		%	NOTE: This method supports array inputs OBJ!
			
			% Handle input arguments
			narginchk(1,1);
			
			% Reverse waypoints
			for i = 1:numel(obj)
				% Cast NBR property to sign preserving data type!
				nbrS = double(obj(i).nbr);
				
				obj(i) = Waypoints(...
					+fliplr(obj(i).x),...
					+fliplr(obj(i).y),...
					-fliplr(obj(i).s) + obj(i).s(end),... % in case of non-equally distributed x/y
					-fliplr(obj(i).curv),...
					+fliplr(obj(i).head) + pi,... 
					+fliplr(obj(i).type),...
				-fliplr(nbrS) + nbrS(end) + 1); 
			end%for
			
		end%fcn
		
		function obj = rotate(obj, phi)
		% ROTATE	Rotate WAYPOINTS object.
		%	OBJ = ROTATE(OBJ,PHI) rotates OBJ by an angle PHI in radians
		%	around the origin.
		% 
		%	OBJ = ROTATE(OBJ) applies the default value -OBJ(1).head(1) for
		%	PHI.
		%	
		%	If OBJ is an array of WAYPOINTS objects, PHI is applied to all
		%	array elements.
		%
		%	NOTE: If OBJ is an array, PHI is applied to all elements of
		%	OBJ.
			
			% Handle input arguments
			narginchk(1, 2);
			
			if nargin < 2
				phi = -obj(1).head(1);
			end%if
			
			if numel(phi) ~= 1 || ~isnumeric(phi)
				error(['Method ROTATE requires a numeric input',...
					' argument with one element.']);
			end%if
			
			
			% Rotate
			% rotation matrix in R^2.
			rotMat = [cos(phi) -sin(phi); 
					  sin(phi) +cos(phi)];
			
			for i = 1:numel(obj)
				% much faster than rotMat*[obj.x;obj.y]
				xy_new = [obj(i).x' obj(i).y']*rotMat';

				obj(i) = Waypoints(...
					xy_new(:, 1)', ...
					xy_new(:, 2)', ...
					obj(i).s, ...
					obj(i).curv, ...
					obj(i).head + phi, ...
					obj(i).type, ...
					obj(i).nbr);
			end%for
			
		end%fcn
		
		function obj = setStartIndex(obj,idx)
		% SETSTARTINDEX		Reorder waypoints.
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
				obj.s(idx:end) - obj.s(idx),...
				obj.s(1:idx-1) + obj.s(end) - obj.s(idx) + delta];
			
			%%% Reorder waypoints
			% Make sure that property NBR is monotonically increasing and
			% starts at 1
			obj = Waypoints(...
				[obj.x(idx:end),	obj.x(1:idx-1)],...
				[obj.y(idx:end),	obj.y(1:idx-1)],...
				s_new,...
				[obj.curv(idx:end),	obj.curv(1:idx-1)],...
				[obj.head(idx:end),	obj.head(1:idx-1)],...
				[obj.type(idx:end),obj.type(1:idx-1)],...
				[obj.nbr(idx:end),	obj.nbr(1:idx-1) + obj.nbr(end)] - obj.nbr(idx) + 1);
			
		end%fcn
		
		function obj = shiftBy(obj, P)
		% SHIFTBY	Shift waypoints by offset.
		%	OBJ = SHIFTBY(OBJ,P) shifts waypoints so that the initial point
		%	is [OBJ.x(1)+P(1) OBJ.y(1)+P(2)].
		% 
		%	NOTE: If OBJ is an array, P is applied to all array elements!
			
			% Handle input arguments
			narginchk(2,2);
			
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
					obj(i).curv,...
					obj(i).head,...
					obj(i).type,...
					obj(i).nbr);
			end%for
			
		end%fcn
		
		function obj = shiftTo(obj, P)
		% SHIFTTO	Shift waypoints to point.
		%	OBJ = SHIFTTO(OBJ,P) shifts waypoints to new inital point P.
		%	
		%	OBJ = SHIFTTO(OBJ) applies the default value [0 0] for P.
		%	
		%	NOTE: If OBJ is an array of WAYPOINTS, P is applied to the
		%	first array element, the others are shifted to maintain their
		%	positioning with respect to the first WAYPOINTS object!
			
			
			% handle input arguments
			narginchk(1,2);
			
			if nargin < 2
				P = [0 0];
			end%if
			
			if numel(P) ~= 2 || ~isnumeric(P)
				error(['Method SHIFTTO requires a numeric input',...
					' argument with two elements.']);
			end%if
			
			% Shift waypoints and maintain positional relation if OBJ is an
			% array
			P0 = obj(1).initPoint;
			for i = 1:numel(obj)
				Pi = P + obj(i).initPoint - P0;
				obj(i) = Waypoints(...
					obj(i).x - obj(i).x(1) + Pi(1),...
					obj(i).y - obj(i).y(1) + Pi(2),...
					obj(i).s,...
					obj(i).curv,...
					obj(i).head,...
					obj(i).type,...
					obj(i).nbr);
			end%for
			
		end%fcn
		
		function obj = unwrap(obj)
		% UNWRAP	Correct heading angle property PHI.
		%	OBJ = UNWRAP(OBJ) applies UNWRAP(MOD(_,2*pi)) to property PHI.
		%	
		%	NOTE: This method supports array inputs OBJ!
		% 
		%	See also UNWRAP, MOD.
			
			for i = 1:numel(obj)
				obj(i) = Waypoints(...
					obj(i).x, ...
					obj(i).y, ...
					obj(i).s, ...
					obj(i).curv, ...
					unwrap( mod(obj(i).head, 2*pi) ), ...
					obj(i).type, ...
					obj(i).nbr);
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
			fid = Waypoints.openFileWithExt(fn);
			
			switch format
				case 'raw'
					fprintf(fid,...
						'%12s %12s %12s %12s %12s %4s %3s \n',...
						'x [m]','y [m]','s [m]','k [1/m]','phi [rad]','type','nbr');
					fprintf(fid,...
						'%+12.3f %+12.3f %12.3f %12.3f %12.3f %4d %3u \n',...
						[obj.x; obj.y; obj.s; obj.curv; obj.head;...
						double(obj.type); double(obj.nbr)]);
					
				case 'CarMaker'
					write2file_CarMaker(fn, obj, varargin{:});
					
				case 'rompac'
					fprintf(fid,...
						'%s %10s %12s %12s %12s %12s \n',...
						'#', 'x_[m]', 'y_[m]', 'phi_[rad]','SD_right_[m]','SD_left_[m]');
					fprintf(fid,...
						'%+12.3f %+12.3f %12.3f %12.3f %12.3f \n',...
						[obj.x; obj.y; obj.head; varargin{1}']);
					
				otherwise
					error('WAYPOINTS:write2File','Unknown FORMAT specifier');
					
			end%switch
			
			% Close file
			fclose(fid);
			
		end%fcn
		
		function [out,latOff_LAD,angDev_LAD,curvat_LAD,isValid_LAD] = ...
				laneTracking(obj,xyCG_global,yawAngle_global,LAD,mode)
		% LANETRACKING  Calculate lane tracking pose.
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
			
			% ensure row format
			LAD = LAD(:)';
			
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
			logIndx1 = bsxfun(@le, xyLAD_T(1,:)' - deltaX_max/2, obj_T.x);
			logIndx2 = bsxfun(@ge, xyLAD_T(1,:)' + deltaX_max/2, obj_T.x);
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
			[indl,indu] = obj.interpIndexRange(numIndCol,[1,obj.nbrel],m);
			
			% preallocation of for-loop variable
			lanePose_LAD_candidates = zeros(3,length(numIndCol));
			try
				% Inter- bzw. Extrapoliere y-Werte der transf. Solltrajektorie
				for i = 1:length(numIndCol)
					% same result like interp1(..,'spline') but faster
					lanePose_LAD_candidates(:,i) = spline(...
						obj_T.x(indl(i):indu(i)),...
						[obj_T.y(indl(i):indu(i));...
						 obj_T.head(indl(i):indu(i));...
						 obj_T.curv(indl(i):indu(i))],...
						xyLAD_T(1, numIndRow(i)));
				end%for
			catch exception
				plotLaneTracking(obj,xyCG_global,yawAngle_global,LAD,numIndCol,obj_T,xyCG_T);
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
	
	
	%%% User-facing methods (plot related)
	methods
		
		function [h,ax] = curvplot(varargin)
		%CURVPLOT	Plot waypoints, heading and curvature.
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
				ax = [subplot(2, 2, [1,3]); subplot(2, 2, 2); ...
					subplot(2, 2, 4)];
			elseif numel(ax) < 3
				error('If specified, three axes handles are required!');
			end%if
			

			h = gobjects(3, 1);
			h(1) = plot(ax(1), obj, opts{:});
			
			h(2) = plot(ax(2), obj.s, obj.head, opts{:});
			grid(ax(2), 'on');
			ylabel(ax(2), 'Heading [rad]');
			
			h(3) = plot(ax(3), obj.s, obj.curv, opts{:});
			grid(ax(3), 'on');
			xlabel(ax(3), 'Length [m]');
			ylabel(ax(3), 'Curvature [1/m]');
			
			% linkaxes(_, 'x') sets axes XLimMode to manual which can cause
			% data of subsequent calls of this method to be not completely
			% shown. Therefore reset axes limit mode to 'auto'.
			linkaxes(ax(2:3), 'x'); 
			set(ax, 'XlimMode', 'auto');
			
		end%fcn
		
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
		%	See also PLOT, WAYPOINTS/CURVPLOT, WAYPOINTS/PLOTTANGENT,
		%	WAYPOINTS/PLOTDIFF.
		
			% Default plot style if unspecified
			opts = {'LineWidth',1};
			
			% Handle accepted input argument configurations
			if nargin == 1
				% Easy case: plot() was called with just a WAYPOINTS object
				obj = varargin{1};
				[h, ax] = plot_raw([], obj, opts{:});
				
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
				[h, ax] = plot_raw(ax, obj, opts{:});
				
			end%if
			
			% Apply plot styles
			grid(ax, 'on');
			axis(ax, 'equal');
			title(ax, getTitleCellString(obj));
			ylabel(ax, 'y [m]');
			xlabel(ax, 'x [m]');
		
		end%fcn
		
		function h = plotdiff(obj, fh)
		%PLOTDIFF	Plot waypoints with specific appearance.
		%	PLOTDIFF(OBJ) plots the different types of waypoints OBJ using
		%	the according pre-defined type-specific color and marker.
		%	
		%	PLOTDIFF(OBJ,FH) marker symbols are switched periodically at
		%	indices specified by the function handle FH. Some usefull
		%	values of FH might be
		%	  (D) @(OBJ) diff(OBJ.type) 
		%	  (.) @(OBJ) diff(OBJ.type) | diff(sign(OBJ.curv))
		%	  (.) @(OBJ) diff(OBJ.nbr)
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
				fh = @(arg) diff(arg.type); 
			end%if
			
			if ~isa(fh, 'function_handle')
				error('Second input argument must be of class function handle!')
			end%if
			
			
			% index numbering of number of elements of OBJ
			indBase = 1:obj.nbrel;
			
			
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
				h(i) = plot_raw([], sd,...
					'LineStyle','none',...
					'Color',	colorStyles{sd.type(1)+2},...
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
			
			n = obj.nbrel;
			xx = (0:n-1)/n;
			xLblString = 'index [1/index_{max}]';
			
			switch prop
				case 'curv'
					yy = obj.curv;
					yLblString = 'curvature [1/m]';
				
				case 'head'
					yy = obj.head;
					yLblString = 'heading [rad]';
					
				case 's'
					yy = obj.s;
					yLblString = 's [m]';
					
				case 'type'
					yy = obj.type;
					yLblString = 'type [-]';
					
				case 'nbr'
					yy = obj.nbr;
					yLblString = 'nbr [-]';
					
				otherwise
					error(['Unknown property string. ',...
						'Type help plotdiff_ for valid property strings.']);
					
			end%switch
			
			% create WAYPOINTS object for plotting
			obj_new = Waypoints(xx, yy, obj.s, obj.curv, obj.head,	obj.type, obj.nbr);
			
			% plot
			h = plotdiff(obj_new, fh);
			
			% apply plot styles
			axis normal
			grid on
			xlabel(xLblString);
			ylabel(yLblString);
			title('');
			
		end%fcn
		
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
			if any(ind > obj.nbrel)
				warning('WAYPOINTS:plottangent:indxOutOfRange', ...
					['Some indexes exceed the number N=%i of waypoints ', ...
					'These tangents will not be plotted.'], obj.nbrel)
				ind(ind > obj.nbrel) = [];
			end%if
			
			% Plot the tangents and the corresponding waypoints
			hold on
			for i = 1:length(ind)
				iind = ind(i);
				
				% Marker of tangent point
				h(i+1,1) = plot(obj.x(iind), obj.y(iind), opts_marker{:}, opts{:});
				
				% Length of tangent
				[r1,r2] = obj.scaleTangentToAxis(xLimits, yLimits, ...
					[obj.x(iind) obj.y(iind)], obj.head(iind));
				
				% Start/end point of tangent
				Pstart = [...
					obj.x(iind)+r2*cos(obj.head(iind)); ...
					obj.y(iind)+r2*sin(obj.head(iind))];
				Pstop = [...
					obj.x(iind)+r1*cos(obj.head(iind)); ...
					obj.y(iind)+r1*sin(obj.head(iind))];
				
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
						sd_T.curv,...
						sd_T.head,...
						sd_T.type,...
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
				h_V = plottangent(sd_T,indx,'k');
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
		
	end%methods
	
	
	methods (Access = private)
		
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
				legend('show', 'Location','best');
			end%if
			
			% Reset axes to initial state
			set(axh, 'NextPlot',npState);
			
		end%fcn
		
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
		
		function PROP = resample_on_s(obj, field, S)
			
			% find indices where property FIELD of OBJ changes
			indPropChange = find( diff(obj.(field)) ~= 0) + 1;
			
			% get corresponding length values
			sPropChange = obj.s(indPropChange);
			
			% initialize resampled data with initial value of original data
			PROP = obj.(field)(1) * ones(size(S), class(obj.(field)));
			
			% set following segments defined by INDPROPCHANGE/SPROPCHANGE
			for i = 1:numel(indPropChange)
				PROP(S > sPropChange(i)) = obj.(field)(indPropChange(i));
			end%for
			
		end%fcn
		
		function cellStr = getTitleCellString(obj)
			
			if numel(obj) > 1
				cellStr = {sprintf('Comparing %d paths', numel(obj))};
			else
				
				if all(obj.type(1) == obj.type)
					typeStr = obj.curvTypes{obj.type(1) + 2};
					typeStr(1) = upper(typeStr(1));
				else
					typeStr = 'Undefined Waypoints';
				end%if
				cellStr	= {...
					sprintf('%s (N=%u;  s=%.2f)', ...
					typeStr, obj.nbrel, obj.s(end)),...
					sprintf(' (%.2g;%.2g) \\rightarrow (%.2g;%.2g)', ...
					obj.x(1), obj.y(1), obj.x(end), obj.y(end)),...
					};
			
			end%if
			
		end%fcn
		
		function write2file_CarMaker(fn, obj, varargin)
			
			% open file
			[fid, fileExt] = Waypoints.openFileWithExt(fn);
			
			switch fileExt
				case '.csv'
					% x .. x-coordinate [m]
					% y .. y-coordinate [m]
					% z .. altitude [m]
					% q .. slope
					% wl/wr .. track width left/right [m]
					% ml/mr .. margin width left/right [m]
					N = obj.nbrel;
					if nargin < 3
						lw_left		= ones(N, 0);
						lw_right	= ones(N, 0);
					else
						lw_left		= varargin{1}(1)*ones(N, 1);
						lw_right	= varargin{1}(2)*ones(N, 1);
					end
					if nargin < 4
						mw_left		= ones(N, 0);
						mw_right	= ones(N, 0);
					else
						mw_left		= varargin{2}(1)*ones(N, 1);
						mw_right	= varargin{2}(2)*ones(N, 1);
					end
					fprintf(fid,':	x	y	z	q	wl	wr	ml	mr\n');
					fprintf(fid,'#\n# IPG ROADDATA\n');
					fprintf(fid,'# This file was created automatically:\n');
					fprintf(fid,'#  Export source: class WAYPOINTS\n');
					fprintf(fid,'#  Export date: %s\n',datestr(now));
					fprintf(fid,'#\n# Add your comments here!\n#\n#\n');
					fprintf(fid,'#	x	y	z	q	wl	wr	ml	mr\n');
					dlmwrite(fn,...
						[obj.x(:), obj.y(:), zeros(N,1), zeros(N,1), ...
						lw_left, lw_right, mw_left, mw_right],...
						'-append', ...
						'delimiter','	', ...
						'precision', '%+f')
						
				case '.kml'
					if isempty(varargin)
						error('WAYPOINTS:WRITE2FILE:UTM_ZoneMissing', ...
							'You need to specify the UTM zone!');
					end%if
					[lat,lon] = utm2ll(obj.x, obj.y, varargin{1});
					
					fprintf(fid, ...
						'<Placemark>\n\t<LineString>\n\t\t<coordinates>\n');
					fprintf(fid, ...
						'\t\t\t%+.12f,%+.12f,%+.12f \n', ...
						[lon(:)'; lat(:)'; zeros(size(lon(:)'))]);
					fprintf(fid, ...
						'\t\t</coordinates>\n\t</LineString>\n</Placemark>');
					
				otherwise
					error('CarMaker export only supports .csv and .kml files!.')
			end%switch
			
		end%fcn
		
	end%methods
	
	
	%%% GET-Methods
	methods
		
		function pos = get.initPoint(obj)
			pos = [obj.x(1) obj.y(1)];
		end%fcn
		
		function pos = get.termPoint(obj)
			pos = [obj.x(end) obj.y(end)];
		end%fcn
		
		function n = get.nbrel(obj)
			n = numel(obj.x);
		end%fcn
		
	end%methods
	
	
	%%% SET-Methods
	methods
		
		function obj = set.s(obj, value)
			% The length can not decrease from one to following point
			if any(diff(value) < 0)
				error('Property S must be monotonically increasing!')
			end%if
			
			obj.s = value;
		end%fcn
		
		function obj = set.type(obj, value)
			
			% Check value range
			if any(value < -1, 2) || any(value > 4, 2)
				error(['Property TYPE out of range! ',...
					'Type help WAYPOINTS/TYPE for valid values.']);
			end%if
			
			% Ensure class
			if ~isa(value, 'int8')
				obj.type = int8(value);
			else
				obj.type = value;
			end%if
			
		end%fcn
		
		function obj = set.nbr(obj, value)
			
			% NBR must be monotonically increasing
			if any(diff(double(value)) < 0)
				error('Input NBR must be monotonically increasing!')
			end%if
			
			% Make sure NBR starts with 1
			if value(1) ~= 1
				fprintf('Warning: Modifying property NBR to start at 1!');
				value = value - value(1) + 1;
			end%if
			
			% Ensure class
			if ~isa(value, 'uint16')
				obj.nbr = uint16(value);
			else
				obj.nbr = value;
			end%if
			
		end%fcn
		
	end%SET-Methods
	
	
	%%% Static-Methods
	methods (Static)
		
		function [r1,r2] = scaleTangentToAxis(xLimits,yLimits,xy,phi)
		%SCALETANGENTTOAXIS		Scale length of tangent to axis limits.
		%   [R1,R2] = SCALETANGENTTOAXIS(XLIMITS,YLIMITS,XY,PHI) calculates
		%   the lengths R1 and R2 for the tangent at point XY with angle
		%   PHI, so that the tangent does not exceed given limits XLIMITS =
		%   [xMin xMax] and YLIMITS = [yMin yMax].
		%	Length R1 counts in the direction of PHI, whereas R2 counts in
		%	the opposite direction.
		%
		%	The intended usage is for plotting something like tangents in
		%	existing plot figures maintaining the current axis limits.
			
			
			% assign inputs to meaningful variables
			xMin = xLimits(1);
			xMax = xLimits(2);
			yMin = yLimits(1);
			yMax = yLimits(2);
			xT = xy(1);
			yT = xy(2);
			
			
			%%% map phi to [-pi,+pi) (plus/minus pi)
			if phi >= pi
				phi_pmPi = phi - 2*pi;
			else
				phi_pmPi = phi;
			end%if
			
			
			%%% reduce problem to angles of [-pi/2,+pi/2]
			if phi_pmPi <= pi/2 && phi_pmPi >= -pi/2
				% keep the original value
				phi_pmPi_ = phi_pmPi;
				isPhiReversed = false;
			else
				% consider the opposite direction of phi which is whithin
				% [-pi/2,+pi/2]
				phi_pmPi_ = mod(phi_pmPi - pi,pi); % ??? consider the sign of x in y for mod(x,y)
				isPhiReversed = true;
			end%if
			
			
			%%% calculate the equation of a line at [xT yT] for phi
			% y(xT) = yT -> k*xT + d = yT -> d = yT - k*xT
			k1 = tan(phi_pmPi_); % tangents slope
			d = yT - k1*xT; % tangents offset
			y_x = @(x) k1*x + d; % tangents equation of a line
			
			
			%%% calculate the plot-range [xL,xR] depending on the sign of
			%%% the tangents slope
			if k1 >= 0
				
				if y_x(xMax) <= yMax
					xR = xMax;
				else
					xR = (yMax-d)/k1;
				end%if
				
				if y_x(xMin) >= yMin
					xL = xMin;
				else
					xL = (yMin-d)/k1;
				end%if
				
			else
				
				if y_x(xMax) >= yMin
					xR = xMax;
				else
					xR = (yMin-d)/k1;
				end%if
				
				if y_x(xMin) <= yMax
					xL = xMin;
				else
					xL = (yMax-d)/k1;
				end%if
				
			end%if
			
			
			%%% calculate the tangent lengths
			r1 = +sqrt((xR-xT)^2 + (y_x(xR)-y_x(xT))^2);
			r2 = -sqrt((xL-xT)^2 + (y_x(xL)-y_x(xT))^2);
			
			% fix tangent lengths if the slope is +-inf
			if phi_pmPi_ == +pi/2
				r1 = yMax - yT;
				r2 = yMin - yT;
			elseif phi_pmPi == -pi/2
				r1 = -(yMin - yT);
				r2 = -(yMax - yT);
			end%if
			
			
			%%% switch R1/R2 and their signs if the opposite problem was
			%%% considered
			if isPhiReversed
				dummy = r1;
				r1 = -r2;
				r2 = -dummy;
			end%if
			
		end%fcn
		
		function [indl,indu,ind] = interpIndexRange(ind,indMinMax,m)
		%INTERPINDEXRANGE	Index range for interpolation.
		%	[INDL,INDU] = INTERPINDEXRANGE(IND,INDMINMAX,M) returns the
		%	lower and upper indices INDL and INDU within index boundaries
		%	[INDMINMAX(1) INDMINMAX(2)] using an index difference M for
		%	given indices IND.
		%	
		%	Indices INDL/INDU define the range of interpolation, basically
		%	IND-M and IND+M but ensure INDL > INDMINMAX(1) and INDU <
		%	INDMINMAX(2).
		
		
			%%% handle input arguments
			% input IND
			if  ~isvector(ind) % check size
				error('Size of IND must be vector!');
			else
				ind = ind(:); % ensure column orientation
			end%if
			
			% input INDMINMAX
			if any(size(indMinMax) ~= [1 2]) % check size
				error('Size of INDMINMAX must be 1-by-2!');
			else
				indMin = indMinMax(1);
				indMax = indMinMax(2);
			end%if
			
			% input M
			if nargin < 3 % apply default value if undefined
				m = 1;
			elseif ~isscalar(m) % check isze
				error('Size of M must be scalar!');
			elseif m < 1 % check value
				error('Value of m must be > 0!');
			end%if
			
			
			%%% check input arguments plausibility
			if any(ind<indMin | ind>indMax)
				error('laneTracking:interpIndexRange:IndexOutOfBounds',...
					'One or more indexes out of range [%i,%i].',indMin,indMax);
			end%if
			if indMax-indMin < 2*m
				error('laneTracking:interpIndexRange:IntervalExtOutOfRange',...
					'Index difference M=%i does not fit into given range [%i,%i].',...
					2*m,indMin,indMax);
			end%if
			
			
			%%% calculate interpolating indexes
			% lower/upper interpolating-index unbounded
			indl = ind-m;
			indu = ind+m;
			
			% bounded to INDMINMAX
			indLU = bsxfun(@plus,[indl,indu],...
				max(zeros(size(ind)),indMinMax(:,1)-indl) + ...
				min(zeros(size(ind)),indMinMax(:,2)-indu));
			indl = indLU(:,1);
			indu = indLU(:,2);
			
		end%fcn
		
		function [fid,fileExt] = openFileWithExt(fn, fExt_set)
			
			[~,fileName,fileExt] = fileparts(fn);
			
			if nargin > 1
				if isempty(fileExt)
					warning('WAYPOINTS:write2file:fileExtension',...
						'Adding file extension ''%s'' to file name!',...
						fExt_set);
				elseif ~strcmp(fileExt, fExt_set)
					warning('WAYPOINTS:write2file:fileExtension',...
						'Replacing file extension ''%s'' by ''%s''!',...
						fileExt, fExt_set);
				end%if
				fn = [fileName, fExt_set];
			end%if
			
			% open file with write-permission
			fid = fopen(fn,'w');
				
		end%fcn
		
		function obj = ll2Waypoints(lat, lon, varargin)%#codegen
		%LL2WAYPOINTS	Convert LAT/LON coordinates to WAYPOINTS.
		%	OBJ = LL2WAYPOINTS(LAT,LON)
			
			% LAT/LON are required, accept two optional name-value pair
			% arguments
			narginchk(2, 6);
			
			% Parse name-value pair arguments
			p = inputParser;
			p.FunctionName = 'll2Waypoints';
			addParameter(p, 'head', [], @isnumeric);
			addParameter(p, 'curv', [], @isnumeric);
			parse(p, varargin{:});
			
			% Convert from lat/lon to UTM
			[x,y] = ll2utm(lat(:)', lon(:)');
			s = sFrom_x_y(x, y);
			
			if isempty(p.Results.head)
				head = curvHeadFrom_x_y(x, y);
			else
				head = p.Results.head(:)';
			end%if
			
			if isempty(p.Results.curv)
				curv = curvFrom_head_s(head, s);
			else
				curv = p.Results.curv(:)';
			end%if
			
			% Create WAYPOINTS object
			obj = Waypoints(x, y, s, curv, head);
			
		end%fcn
		
		function obj = pp2Waypoints(t, ppx, ppy, mode)
		% PP2WAYPOINTS	Convert piecewise polynomial structure to WAYPOINTS. 
		%   OBJ = PP2WAYPOINTS(T,PPX,PPY) creates the instace OBJ of class
		%   WAYPOINTS from piecewise polynomial structures PPX and PPY sampled
		%   at T. 
		% 
		%	See also SPLINE.
			
			narginchk(2, 4);
			
			if nargin < 4
				mode = 'int';
			end%if
			
			% Get dimension via UNMKPP to check for valid piecewise
			% polynomial struct.
			[~,~,~,~,dimX] = unmkpp(ppx);
			
			% Make sure T is a row vector. As a result, PPVAL will return a
			% matrix whose number of columns matches the number of elements
			% of T.
			if ~isvector(t) || ~isrow(t)
				error('Query point must be a row vector!');
			end%if
			
			
			if (nargin < 3) || isempty(ppy)
				% In the case of a two-argument call, PPX needs to return a
				% 2-D array (x and y).
				if dimX ~= 2
					error('Piecewise polynomial struct must have a dimensionality of 2!');
				end%if
				
				% Evaluate piecewise polynomial structs
				xy	 = ppval(ppx, t);
				d1xy = ppval(ppdiff(ppx,1), t);
				d2xy = ppval(ppdiff(ppx,2), t);
				
				% Function handle for path length integration
				ppdx = ppdiff(ppx);
				fun = @(t) sqrt(sum(ppval(ppdx, t).^2, 1));
				
			else
				% In the case of a three-argument call, both piecewise
				% polynomial structs must return a 1-D array.
				[~,~,~,~,dimY] = unmkpp(ppy);
				if (dimX > 1) || (dimY > 1)
					error('Piecewise polynomial structs must have a dimensionality of 1!');
				end%if
				
				% Evaluate piecewise polynomial structs
				xy	 = [ppval(ppx, t); ppval(ppy, t)];
				d1xy = [ppval(ppdiff(ppx,1), t); ppval(ppdiff(ppy,1), t)];
				d2xy = [ppval(ppdiff(ppx,2), t); ppval(ppdiff(ppy,2), t)];
				
				% Function handle for path length integration
				ppdx = ppdiff(ppx);
				ppdy = ppdiff(ppy);
				fun = @(t) sqrt(ppval(ppdx, t).^2 + ppval(ppdy, t).^2);
				
			end%if
			
			switch mode
				case 'cumsum'
					s = [0, cumsum(sqrt( sum(diff(xy,1,2).^2, 1) ))];
					
				otherwise
					% Numerical integration of path length
					s = zeros(1, numel(t));
					for i = 1:numel(t)-1
						s(i+1) = integral(@(t) fun(t), t(i), t(i+1));
					end%for
					s = cumsum(s);
			end%%switch
			
			head = unwrap(atan2(d1xy(2,:), d1xy(1,:)));
			curv = (d1xy(1,:).*d2xy(2,:) - d2xy(1,:).*d1xy(2,:)) ./ ...
					sum(d1xy.^2, 1).^(3/2);
			
			obj = Waypoints(xy(1,:), xy(2,:), s, curv, head);
			
		end%fcn
		
		function obj = xy2Waypoints(x, y)
		%XY2WAYPOINTS	Convert x/y coordinates to WAYPOINTS. 
		%   OBJ = XY2WAYPOINTS(X,Y) creates the instace OBJ of class WAYPOINTS
		%   from X/Y-coordinates. X and Y must be vectors with the same
		%   number of elements!
		% 
		%   OBJ = XY2WAYPOINTS(XY) X/Y-coordinates are given in the array XY
		%   of size 2-by-n.
		%
		%	Property S is calculated using cumulative sum.
		%	Properties PHI and CURV are calculated as difference quotients.
		%	Property TYPE is set to -1.
		%	Property NBR is set to 1.
		
		
			narginchk(1, 2);
		
			%%% Handle input arguments
			if nargin < 2
				if size(x, 1) ~= 2
					error('When using one argument syntax, array must have 2 rows!');
				else
					y = x(2, :);
					x = x(1, :);
				end%if
			end%if
			
			if ~isvector(x) || ~isvector(y)
				error('Input arguments X/Y must be vectors!');
			end%if
			
			if numel(x) ~= numel(y)
				error('Input arguments X/Y must have the same number of elements!');
			else
				x = x(:)';
				y = y(:)';
			end%if
			
			
			% Check for NaN's: GRADIENT will fail to return a numeric value
			% for terminal points in case of NANs
			isNan = or(isnan(x), isnan(y));
			
			%%% Convert to WAYPOINTS
			s = sFrom_x_y(x(~isNan), y(~isNan));
			[head,~] = curvHeadFrom_x_y(x(~isNan), y(~isNan));
			curv = curvFrom_head_s(head, s);
			
			obj = Waypoints(x(~isNan), y(~isNan), s, curv, head, -1, 1);
			
		end%fcn

	end%methods
	
end%class
