classdef segDat
%SEGDAT 	Street segment data.
%	
%	OBJ = SEGDAT(X,Y,S,K,PHI,TYPE,NBR) creates the street segment object
%	OBJ of class SEGDAT representing the geometric data of a street segment
%	using the properties listed below.
%	
%	A curvature K > 0 (K < 0) indicates a left (right) turn by moving from
%	[X(1),Y(1)] -> [X(end),Y(end)].
%	
%	NBR is usually 1, but becomes interesting when street segments are
%	connected together. Then NBR represents the number of the street
%	segment in the path of connected segments.
%	
%	SEGDAT Properties:
%	 X		- x-coordinate [m].
%	 Y		- y-coordinate [m].
%	 S		- Length measured from the starting point [x(1);y(1)] [m].
%	 K		- Curvature [1/m].
%	 PHI	- Angle beetween the positive x-axis and the tangent [rad].
%	 TYPE	- Integer indicating the street segment type.
%	 NBR	- Street segment number of distinctive segments.
%	
%	SEGDAT Methods:
%	 - DESIGN
%	 changeSignOfCurvature - Change street segments curvature sign.
%	 plus	- Connect street segments using '+'.
%	 reverseDirection - Reverse street segment direction.
%	 rotate	- Rotate street segment.
%	 selectIndexRange - Select subset of street segment.
%	 setStartIndex - Set an new starting starting element.
%	 shiftBy - Shift street segment by point.
%	 shiftTo - Shift street segment to point.
%	 
%	 - ANALYSIS
%	 plot		 - Plot street segments.
%	 plotcmp	 - Plot variable number of street segments to same figure.
%	 plotdiff	 - Plot street segments with specific appearance.
%	 plotdiff_	 - Plot street segment property with specific appearance.
%	 plottangent - Plot street segments and specified tangents.
%	 
%	 - MISC
%	 laneTracking - Get the lane tracking pose.
%	 write2file	- Write segment data to file.
%	 xy2segDat - Convert x/y-coordinates to SEGDAT object.
%	 
%	
%	See also LKASEGMENT, LKASEGMENTSTRAIGHT, LKASEGMENTCIRCLE,
%	LKASEGMENTCLOTHOID.
%

% DEVELOPMENT NOTES:
%	(1) Add methods for plot labeling (x/y labels, title, ...). DONE
%	(2) Show the clothoid direction in plots.

% Subject: lka
% Author: $Author$
% Date: $LastChangedDate$
% Revision: $Revision$



	properties (SetAccess = private)
		% 1-by-n double of X-coordinate [m].
		x
		
		% 1-by-n double of Y-coordinate [m].
		y
		
		% 1-by-n double of Covered distance [m].
		s
		
		% 1-by-n double of Curvature [1/m].
		k
		
		% 1-by-n double of Tangent angle [rad].
		phi
		
		% 1-by-n int8 of curvature type [-]:
		%	-1 .. unknown
		%	 0 .. straight
		%	 1 .. circular
		%	 2 .. clothoid
		type
		% See also SEGDAT/CURVTYPES
		
		% 1-by-n uint16 of Street segment number [-].
		nbr
	end%properties
	
	
	properties (Constant, Hidden)
		
		% curvature types
		curvTypes = {'unknown','straight','circle','clothoid'};
		
		% Plot color: one color per SEGDAT/TYPE
		plotColor = {'m','b','g','r'};
		
		% Plot marker symbol/size
		plotMarker = {...
			'o','diamond';... % marker symbol
			5,	4}; % marker size
		
	end%properties
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%% CONSTRUCTOR
	methods
		
		function obj = segDat(x,y,s,k,phi,type,nbr)
		%SEGDAT 	Create an instance of the SEGDAT object.
		
		
			% all inputs have to be non-empty column vectors
			[m,n] = size(x);
			
			% expand the lenght of TYPE/NBR if they are scalar
			if isscalar(type)
				type = type*ones(m,n);
			end%if
			if isscalar(nbr)
				nbr = nbr*ones(m,n);
			end%if
			
			% just check X and check the size of all others against X
			if (m>1) || (n<1)
				error('SEGDAT:segDat:nonemptyColumnVectors',...
					'Inputs have to be non-empty row vectors')
			end%if
			if ~isequal(size(x),size(y),size(s),size(k),size(phi),size(type),size(nbr))
				error('SEGDAT:segDat:unequalInputArgumentSizes',...
					'Inputs must have the same size.')
			end%if
			
			% check class of TYPE
			if ~isa(type,'int8')
				type = int8(type);
			end%if
			
			% check class of NBR
			if ~isa(nbr,'uint16')
				nbr = uint16(nbr);
			end%if
			
			obj.x		= x;
			obj.y		= y;
			obj.s		= s;
			obj.k		= k;
			obj.phi		= phi;
			obj.type	= type;
			obj.nbr		= nbr;
			
		end%fcn
		
	end%CONSTRUCTOR-methods
	
	
	%%% User-facing methods
	methods
		
		function obj = changeSignOfCurvature(obj)
		% CHANGESIGNOFCURVATURE		Change street segments curvature sign.
		%	OBJ = CHANGESIGNOFCURVATURE(OBJ) changes the curvature sign of
		%	street segment OBJ while maintaining the initial point.
		%	
		%	The other properties are manipulated accordingly.
		
			% get starting point/angle
			P0 = [obj.x(1); obj.y(1)];
			phi0 = obj.phi(1);
			
			% shift to origin and rotate so initial slope is zero
			obj = shiftTo(obj);
			obj = rotate(obj,-phi0);
			
			% invert the curvature
			obj = segDat(...
				+obj.x,...
				-obj.y,...
				+obj.s,...
				-obj.k,...
				-obj.phi,...
				+obj.type,...
				+obj.nbr);
			
			% undo the shift/rotate procedure
			obj = rotate(obj,phi0);
			obj = shiftTo(obj,P0);
			
		end%fcn
				
		
		function [lkaSeg,e,c,d] = fitStraight(obj,indMinMax,doPlot)
		%FITSTRAIGHT	Fit a straight line to a set of points.
		%	[LKASEG,E,C,D] = FITSTRAIGHT(OBJ) fits a straight
		%	line with slope C and initial offset D to a set of given points
		%	(OBJ.X,OBJ.Y) minimizing the error E. Object LKASEG is of class
		%	LKASEGMENTSTRAIGHT
		%	
		%	[LKASEG,E,C,D] = FITSTRAIGHT(OBJ,INDMINMAX) lets you specify a
		%	start index I = INDMINMAX(1) and end index U = INDMINMAX(2) for
		%	considering only the range I:U for the fitting procedure.
		%	
		%	[LKASEG,E,C,D] = FITSTRAIGHT(OBJ,INDMINMAX,DOPLOT) allows to
		%	disable the plot for checking the fitting result visually,
		%	which is enabled by default.
		%	
		%	Parameters C and D model the fitted line according to
		%	  y_fit = C*OBJ.X + D
		%	
		%	Minimization is performed in the least-squares sense minimizing
		%	the sum of squared errors:
		%	  SUM[(Y(i)-C*X(i) - D)^2]
		%	   i
		%	
		%	See also LKASEGMENTSTRAIGHT.
			
			%%% handle input arguments
			if nargin < 2
				indMin = 1;
				indMax = numel(obj.x);
			else
				indMin = indMinMax(1);
				indMax = indMinMax(2);
			end%if
			
			if nargin < 3
				doPlot = true;
			end%if
			
			% extract relevant x/y data
			xsub = obj.x(indMin:indMax);
			ysub = obj.y(indMin:indMax);
			
			% create (overdetermined) system of equations
			A = [xsub',ones(indMax-indMin+1,1)];
			b =  ysub';
			
			% solve system of equations: A*[c;d] = b, where y = c*x+d
			cd = (A'*A)\A'*b;
			c = cd(1);
			d = cd(2);
			
			% the fitted y coordinates
			yfit = c*xsub + d;
			
			% create LKASEGMENTSTRAIGHT object
			dx = xsub(end) - xsub(1);
			dy = yfit(end) - yfit(1);
			lkaSeg = lkaSegmentStraight(...
				min(diff(sqrt(xsub.^2 + ysub.^2))),...
				sqrt(dx^2 + dy^2),...
				atan2(dy,dx));
			lkaSeg = shift(lkaSeg,[xsub(1);ysub(1)]);
			
			% calculate error
			e = 1/numel(ysub)*sum((ysub - yfit).^2);
			
			% plot if required
			if doPlot
				plot(obj,'r');
				hold on
				plot(lkaSeg,'b');
			end%if
			
			
		end%fcn
		
		
		function [lkaSeg,e,xc,yc,R] = fitCircle(obj,indMinMax,doPlot)
		%FITCIRCLE	Fit a circle to a set of points.
		%	[LKASEG,E,XC,YC,R] = FITCIRCLE(OBJ) fits a circle
		%	with center(XC,YC) and radius R to a set of given points
		%	(OBJ.X,OBJ.Y) minimizing the error E. Object LKASEG is of class
		%	LKASEGMENTCIRCLE.
		%	
		%	[LKASEG,E,XC,YC,R] = FITCIRCLE(OBJ,INDMINMAX) lets you specify
		%	a start index I = INDMINMAX(1) and end index U = INDMINMAX(2)
		%	for considering only the range I:U for the fitting procedure.
		%	
		%	[LKASEG,E,XC,YC,R] = FITCIRCLE(OBJ,INDMINMAX,DOPLOT) allows to
		%	disable the plot for checking the fitting result visually,
		%	which is enabled by default.
		%	
		%	Minimization is performed in the least-squares sense minimizing
		%	the sum of squared errors:
		%	  SUM[(R(i)^2-R^2)^2]
		%	   i
		%	
		%	See also LKASEGMENTCIRCLE.
			
			%%% handle input arguments
			if nargin < 2 || isempty(indMinMax)
				indMin = 1;
				indMax = numel(obj.x);
			else
				indMin = indMinMax(1);
				indMax = indMinMax(2);
			end%if
			
			if indMin >= indMax
				error('SEGDAT:FITCIRCLE:index',...
					'Start index must be smaller than end index!');
			end%if
			
			if nargin < 3
				doPlot = true;
			end%if
			
			
			% extract relevant x/y data
			xsub = obj.x(indMin:indMax);
			ysub = obj.y(indMin:indMax);
			
			method = 'Kasa';
			switch method
				case 'Kasa'
					[xc,yc,R,e] = fitCircle_Kasa(xsub,ysub);
				otherwise
					% 
			end%switch
			
			% create LKASEGMENTCIRCLE object
			lkaSeg = lkaSegmentCircle([],0,2*pi,R);
			lkaSeg = shift(lkaSeg,[xc+lkaSeg.radius;yc]);
			
			
			% plot if required
			if doPlot
				plot(obj,'r','Marker','o','DisplayName','Original path');
				hold on
				plot(xsub,ysub,'r','Marker','.','DisplayName','Path under test');
				plot(lkaSeg,'b','DisplayName','Fitted path');
			end%if			
			
		end%fcn
		
		
		function obj = plus(obj1,obj2)
		%+ Plus.
		%	OBJ12 = OBJ1 + OBJ2 adds the street segment data OBJ2 with its
		%	starting point to the end point of street segment data OBJ1
		%	resulting in the street segment data OBJ12. 
		%	
		%	End point of OBJ1 and starting point of OBJ2 are assumed to
		%	overlap! Therefore, the resulting segment OBJ consists of one
		%	element less than the total number of elements of OBJ1 and
		%	OBJ2.
		%	
		%	Note that here plus (+) is a non-commutative operation!
			
			obj = segDat(...
				[obj1.x,	obj1.x(end) + obj2.x(2:end) - obj2.x(1)],...
				[obj1.y,	obj1.y(end) + obj2.y(2:end) - obj2.y(1)],...
				[obj1.s,	obj1.s(end) + obj2.s(2:end)],...
				[obj1.k,	obj2.k(2:end)],...
				[obj1.phi,	obj2.phi(2:end)],...
				[obj1.type, obj2.type(2:end)],...
				[obj1.nbr,	obj2.nbr(2:end)+obj1.nbr(end)]);
			
		end%fcn
		
		
		function obj = reverseDirection(obj)
		% REVERSEDIRECTION	Reverse street segment direction.
		%	OBJ = REVERSEDIRECTION(OBJ) reverses the direction of street
		%	segment OBJ so [x(end),y(end)] becomes [x(1),y(1)] and so on.
		%	
		%	The other properties are manipulated accordingly.
			
			%%% handle input arguments
			narginchk(1,1);
			
			%%% reverse direction of segment
			obj = segDat(...
				+flipud(obj.x),...
				+flipud(obj.y),...
				+obj.s(end)-flipud(obj.s),... % in case of non-equally distributed x/y
				-obj.k,...
				+flipud(obj.phi) + pi,... 
				+obj.type,...
				+obj.nbr); 
			
		end%fcn
		
		
		function obj = rotate(obj,phi)
		% ROTATE	Rotate street segment.
		%	OBJ = ROTATE(OBJ,PHI) rotates the street segment OBJ by an
		%	angle PHI in radians.
			
			%%% handle input arguments
			narginchk(2,2);
			
			if numel(phi) ~= 1 || ~isnumeric(phi)
				error(['Method ROTATE requires a numeric input',...
					' argument with one element.']);
			end%if
			
			
			%%% rotate segment
			% rotation matrix in R^2.
			rotMat = [...
				cos(phi) -sin(phi);...
				sin(phi) +cos(phi)];
			
			% much faster than rotMat*[obj.x;obj.y]
			xy_new = [obj.x' obj.y']*rotMat';
			
			obj = segDat(...
				xy_new(:,1)',...
				xy_new(:,2)',...
				obj.s,...
				obj.k,...
				obj.phi+phi,...
				obj.type,...
				obj.nbr); 
			
		end%fcn
		
		
		function obj = selectIndexRange(obj,indRange)
		% SELECTINDEXRANGE	Select subset of street segment.
		%	OBJ = SELECTINDEXRANGE(OBJ,INDRANGE) selects a subset of street
		%	segment OBJ specified by indices INDRANGE.
		%	
		%	Indices INDRANGE must be strictly increasing, otherwise a
		%	warning will be issued!
		%	
		%	If INDRANGE is scalar, it is used as the start index and the
		%	end index is set to the lenght of OBJ.
			
			%%% handle input arguments
			narginchk(2,2);
			
			if isempty(indRange)
				error('SEGDAT:selectIndexRange',...
					['Input argument INDRANGE must not be empty! ',...
					'Type help segDat/selectIndexRange.']);
			end%if
			
			if numel(indRange) < 2
				ind = indRange:numel(obj.x);
			else
				ind = indRange;
			end%if
			
			if any(diff(ind) <= 0)
				warning('segDat:selectIndexRange',...
					['Indices INDRANGE are not strictly increasing! ',...
					'This might result in incorrect results (property S)!']);
			end%if
			
			
			%%% get subset of street segment
			obj = segDat(...
				obj.x(ind),...
				obj.y(ind),...
				obj.s(ind) - obj.s(ind(1)),...
				obj.k(ind),...
				obj.phi(ind),...
				obj.type(ind),...
				obj.nbr(ind));
			
		end%fcn
		
		
		function obj = setStartIndex(obj,indx)
		% SETSTARTINDEX		Reorder street segment.
		%	OBJ = SETSTARTINDEX(OBJ,INDX) reorders the street segment OBJ
		%	so that index INDX becomes the first element.
		%
		%	WARNING: this method only makes sense for circuits!
			
			%%% handle input arguments
			if indx == 1
				return;
			end%
			
			% distance from last to first element
			delta = sqrt((obj.x(1)-obj.x(end))^2 + (obj.y(1)-obj.y(end))^2);
			
			% the length should be stricly monotonically increasing
			s_new = [...
				obj.s(indx:end) - obj.s(indx),...
				obj.s(1:indx-1) + obj.s(end) - obj.s(indx) + delta];
			
			%%% reorder street segment
			obj = segDat(...
				[obj.x(indx:end),	obj.x(1:indx-1)],...
				[obj.y(indx:end),	obj.y(1:indx-1)],...
				s_new,...
				[obj.k(indx:end),	obj.k(1:indx-1)],...
				[obj.phi(indx:end),	obj.phi(1:indx-1)],...
				[obj.type(indx:end),obj.type(1:indx-1)],...
				[obj.nbr(indx:end),	obj.nbr(1:indx-1)]);
			
		end%fcn
		
		
		function obj = shiftBy(obj,P)
		% SHIFTBY	Shift street segment by given point.
		%	OBJ = SHIFTBY(OBJ,P) shifts the street segment OBJ so that its
		%	starting point is [OBJ.x(1)+P(1) OBJ.y(1)+P(2)].
			
			%%% handle input arguments
			narginchk(2,2);
			
			if numel(P) ~= 2 || ~isnumeric(P)
				error(['Method SHIFT requires a numeric input',...
					' argument with two elements.']);
			end%if
			
			
			%%% shift segment
			obj = segDat(...
				obj.x + P(1),...
				obj.y + P(2),...
				obj.s,...
				obj.k,...
				obj.phi,...
				obj.type,...
				obj.nbr); 
			
		end%fcn
		
		
		function obj = shiftTo(obj,P)
		% SHIFTTO	Shift street segment to desired point.
		%	OBJ = SHIFTTO(OBJ,P) shifts the street segment OBJ so that its
		%	starting point [OBJ.x(1) OBJ.y(1)] matches P.
		%	
		%	OBJ = SHIFTTO(OBJ) applies the default value [0 0] for P.
			
			
			%%% handle input arguments
			narginchk(1,2);
			
			if nargin < 2
				P = [0 0];
			end%if
			
			if numel(P) ~= 2 || ~isnumeric(P)
				error(['Method SHIFT requires a numeric input',...
					' argument with two elements.']);
			end%if
			
			
			%%% shift segment
			obj = segDat(...
				obj.x - obj.x(1) + P(1),...
				obj.y - obj.y(1) + P(2),...
				obj.s,...
				obj.k,...
				obj.phi,...
				obj.type,...
				obj.nbr); 
			
		end%fcn
				
		
		function obj = shift(obj,P)
			
			warning('SEGDAT:SHIFT',[...
				'Method SHIFT will be removed in a future release! ',...
				'Use method SHIFTTO instead.']);
			obj = shiftTo(obj,P);
		
		end%fcn
		
		
		function write2file(obj,fn,format,varargin)
		%WRITE2FILE		Write street segment to file.
		%	WRITE2FILE(OBJ,FN) writes the street segment data OBJ to file
		%	with filename FN (also specify extension!).
		%	
		%	WRITE2FILE(OBJ,FN,FORMAT,VARARGIN) lets you specify the output
		%	format. Supported values are:
		%	  - 'raw' (default value)
		%		 One column per property of OBJ. Property names are used as
		%		 labels in the first row.
		%		 
		%	  - 'CarMaker4.0'
		%		 ASCII file (cartesian coordinates), supports at least
		%		 CarMaker v4.0 - v5.0.1.
		%		 Set the lane widths using VARARGIN = [width_left width
		%		 right].
		%		 
		%	  - 'CarMaker4.0_KML'
		%		 KML file (WGS84-coordinates), supports at least
		%		 CarMaker v4.0 - v5.0.1. Use VARARGIN to specify a scalar
		%		 or vector of UTM zone integers.
		%		 The file extension is set to '.kml'.
			
		
			% set the default FORMAT
			if nargin < 3 || isempty(format)
				format = 'raw';
			end%if
			
			
			switch format
				case 'raw'
					% open file
					fid = segDat.openFileWithExt(fn);
					
					fprintf(fid,...
						'%12s %12s %12s %12s %12s %4s %3s \n',...
						'x [m]','y [m]','s [m]','k [1/m]','phi [rad]','type','nbr');
					fprintf(fid,...
						'%+12.3f %+12.3f %12.3f %12.3f %12.3f %4u %3u \n',...
						[obj.x;obj.y;obj.s;obj.k;obj.phi;...
						double(obj.type);double(obj.nbr)]);
					
				case 'CarMaker4.0'
					% open file
					fid = segDat.openFileWithExt(fn);
					
					% x .. x-coordinate [m]
					% y .. y-coordinate [m]
					% z .. altitude [m]
					% q .. slope
					% wl/wr .. track width left/right [m]
					% ml/mr .. margin width left/right [m]
					N = numel(obj.x);
					if nargin < 4
						laneWidth_left	= 0;
						laneWidth_right = 0;
					else
						laneWidth_left	= varargin{1}(1);
						laneWidth_right = varargin{1}(2);
					end
					fprintf(fid,':	x	y	z	q	wl	wr	ml	mr\n');
					fprintf(fid,'#\n# IPG ROADDATA\n');
					fprintf(fid,'# This file was created automatically:\n');
					fprintf(fid,'#  Export source: class SEGDAT\n');
					fprintf(fid,'#  Export date: %s\n',datestr(now));
					fprintf(fid,'#\n# Add your comments here!\n#\n#\n');
					fprintf(fid,'#	x	y	z	q	wl	wr	ml	mr\n');
					dlmwrite(fn,...
						[obj.x(:),obj.y(:),zeros(N,1),zeros(N,1),...
						laneWidth_left*ones(N,1),laneWidth_right*ones(N,1)],...
						'-append',...
						'delimiter','	',...
						'precision','%+f')
					
				case 'CarMaker4.0_KML'
					% open file
					fid = segDat.openFileWithExt(fn,'.kml');
					
					if isempty(varargin)
						error('SEGDAT:WRITE2FILE:UTM_ZoneMissing',...
							'You need to specify the UTM zone!');
					end%if
					[lat,lon] = utm2ll(obj.x,obj.y,varargin{1});
					
					fprintf(fid,...
						'<Placemark>\n\t<LineString>\n\t\t<coordinates>\n');
					fprintf(fid,...
						'\t\t\t%+f,%+f,%+f \n',...
						[lon;lat;zeros(size(obj.x))]);
					fprintf(fid,...
						'\t\t</coordinates>\n\t</LineString>\n</Placemark>');
					
				otherwise
					error('segDat:write2File','Unknown FORMAT specifier');
					
			end%switch
			
			% close file
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
		%	with respect to the street segment OBJ. The logical vector
		%	ISVALID_LAD indicates if the corresponding entry of all other
		%	output arguments is valid (true) or not (false).
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
			
			
			%%% get indices of potential elements of desired path
			% maximaler Abstand zwischen x-Werten der Sollbahn
			deltaX_max = max(abs(diff(obj_T.x)));
			
			% Depending on the shape of the desired path (e.g. closed
			% path), there might be multiple elements of the desired path
			% whose x-coordinates are within the range of LAD +-
			% DELTAX_MAX/2.
			% 
			% Get logical indices where OBJ_T.X is in the range of LADs
			% x-coordinate
			%  rows: LAD 
			%  columns: street segments x-coordinate
			logIndx1 = bsxfun(@le,xyLAD_T(1,:)' - deltaX_max/2,obj_T.x);
			logIndx2 = bsxfun(@ge,xyLAD_T(1,:)' + deltaX_max/2,obj_T.x);
			logIndx = logIndx1 & logIndx2;
			
			% error if no element of LOGINDX is true
			if ~any(any(logIndx))
				plotLaneTracking(obj,xyCG_global,yawAngle_global,LAD,[],obj_T,xyCG_T);
				error('segDat:laneTracking',...
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
			[indl,indu] = obj.interpIndexRange(numIndCol,[1,length(obj.x)],m);
			
			% preallocation of for-loop variable
			lanePose_LAD_candidates = zeros(3,length(numIndCol));
			try
				% Inter- bzw. Extrapoliere y-Werte der transf. Solltrajektorie
				for i = 1:length(numIndCol)
					% same result like interp1(..,'spline') but faster
					lanePose_LAD_candidates(:,i) = spline(...
						obj_T.x(indl(i):indu(i)),...
						[obj_T.y(indl(i):indu(i));...
						 obj_T.phi(indl(i):indu(i));...
						 obj_T.k(indl(i):indu(i))],...
						xyLAD_T(1,numIndRow(i)));
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
		
		function h = plot(obj,varargin)
		%PLOT	Plot the street segment.
		%	PLOT(OBJ) plots OBJ.y over OBJ.x
		%	
		%	PLOT(OBJ,S) additionally applies the line specification
		%	S.
		%	
		%	H = PLOT(...) returns the handle H to lineseries objects.
		%	
		%	The line specification S is a character string supported by the
		%	standard PLOT command. For example
		%		PLOT(OBJ,'LineWidth',2,'Color',[.6 0 0]) 
		%	will create a plot with a dark red line width of 2 points.
		%	
		%	See also PLOT, SEGDAT/PLOTTANGENT, SEGDAT/PLOTDIFF.
		
		
			% apply plot options if unspecified
			if isempty(varargin)
%				varargin = {'o','MarkerSize',2,'MarkerFaceColor','blue'};
				varargin = {'-b','LineWidth',2};
			end%if
			
			% plot street segment
			h = plot_raw(obj,varargin{:});
			
			% apply plot styles
			grid on;
			axis equal;
			title(getLegendCellString(obj));
			ylabel('y [m]');
			xlabel('x [m]');
		
		end%fcn
		
		
		function h = quiver(obj,varargin)
		%QUIVER		Quiver plot of street segment.
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
				varargin = {'-b','LineWidth',2};
			end%if
			
			% plot street segment
			h = quiver_raw(obj,varargin{:});
			
			% apply plot styles
			grid on;
			axis equal;
			title(getLegendCellString(obj));
			ylabel('y [m]');
			xlabel('x [m]');
		
		end%fcn
		
		
		function h = plotcmp(obj,varargin)
		%PLOTCMP	Compare street segments.
		%	PLOTCMP(OBJ,VARARGIN) plots the street segments OBJ and
		%	VARARGIN to the same figure to compare them against each other.
		%	
		%	H = PLOTCMP(...) returns the handle H to lineseries objects.
		%	
		%	See also segDat/plot.
		
			
			h = zeros(1+numel(varargin),1);
			h(1) = plot_raw(obj);
			
			for i = 1:numel(varargin)
				if i == 1; hold all; end
				
				h(i+1) = plot_raw(varargin{i});
				
				if i == numel(varargin); hold off; end
			end%for
			
			% apply plot styles
			grid on;
			axis equal;
			ylabel('y [m]');
			xlabel('x [m]');
		
		end%fcn
		
		
		function h = plotdiff(obj,fh)
		%PLOTDIFF	Plot the street segment with specific appearance.
		%	PLOTDIFF(OBJ) plots each street segment type of OBJ using the
		%	according pre-defined type-specific color and marker.
		%	
		%	PLOTDIFF(OBJ,FH) similar to PLOTDIFF(OBJ) but marker symbols
		%	replace the solid plot line and the marker symbols are switched
		%	periodically at indices specified by the function handle FH.
		%	Some usefull values of FH might be
		%	  (D) @(OBJ) diff(OBJ.type) 
		%	  (.) @(OBJ) diff(OBJ.type) | diff(sign(OBJ.k))
		%	  (.) @(OBJ) diff(OBJ.nbr)
		%	where (D) is used by PLOTDIFF(OBJ). You can also
		%	specify the indices where marker symbols change manually by
		%	using something like
		%		.) @(OBJ) [0 0 1 0 1 0 0 0 1 ...].
		%	Internally FH is evaluated to FH(OBJ) ~= 0.
		%	
		%	H = PLOTDIFF(...) returns a column vector of handles to
		%	lineseries objects, one handle per plotted street segment.
		%	
		%	To modify the plot colors set the hidden property plotColour.
		%	The default colors are blue (straight), green (circle) and red
		%	(clothoid).
		%	
		%	To modify the plot marker set the hidden property plotMarker,
		%	where row one applies to the marker symbol and row two to the
		%	marker size. The default marker symbols are 'o' and 'diamond'
		%	repeated periodically.
		%	
		%	See also SEGDAT/PLOT, SEGDAT/PLOTTANGENT.
		
		
			%%% handle input arguments
			if nargin < 2 || isempty(fh)
				fh = @(arg) diff(arg.type); 
			end%if
			
			if ~isa(fh,'function_handle')
				error('Second input argument must be of class function handle!')
			end%if
			
			
			%%% index numbering of number of elements of OBJ
			indBase = 1:length(obj.x);
			
			
			%%% get the segment-grouping indices and the number of indices
			flag = fh(obj) ~= 0;
			if length(flag) > length(indBase)
				error(['The length of your evaluated function handle ',...
					'FH exceeds the length of the street segment: ',...
					'%u > %u.'],length(flag),length(indBase));
			end%if
			ind = [0,indBase(flag),indBase(end)];
			nbrInd = length(ind);
			
			
			%%% extend the plot settings to the number of segments to plot
			n = ceil(nbrInd/size(obj.plotMarker,2));
			plotMarker_ = repmat(obj.plotMarker,1,n);
			
			
			%%% plot the segments with according plot options
			h = zeros(nbrInd-1,1);
			for i = 1:nbrInd-1
				
				% plot range
				indRange = ind(i)+1:ind(i+1);
				
				% create the segment data to plot
				sd = segDat(...
					obj.x(indRange),...
					obj.y(indRange),...
					obj.s(indRange),...
					obj.k(indRange),...
					obj.phi(indRange),...
					obj.type(indRange),...
					ones(size(indRange))...
					);
				
				if i > 1
					hold on
				end
				h(i) = plot_raw(sd,...
					'LineStyle','none',...
					'Color',obj.plotColor{sd.type(1)+2},...
					'Marker',plotMarker_{1,i},...
					'MarkerSize',plotMarker_{2,i});
				hold off
				
			end%for
			
			% unsure about usefulness
			% line style plotting if no additional ....
% 			if nargin < 2
% 				set(h,'LineStyle','-','LineWidth',2,'Marker','none');
% 			end%if
			
			% apply plot styles
			grid on;
			axis equal;
			title(getLegendCellString(obj));
			ylabel('y [m]');
			xlabel('x [m]');
			
		end%fcn
		
		
		function h = plotdiff_(obj,fh,prop)
		%PLOTDIFF_	Plot street segment property with specific appearance.
		%	PLOTDIFF_(OBJ,FH,PROP) works similar to PLOTDIFF but plots the
		%	property PROP.
		%	
		%	PROP can take the following values:
		%	 'k': plots the curvature over curve length
		%	 'phi': plots the angle over curve length
		%	 's': plots curve length over normed index
		%	
		%	FH can be specified as []. In this case, the default value is
		%	used.
		%	
		%	See also SEGDAT/PLOTDIFF.
			
			n = length(obj.x);
			xx = (0:n-1)/n;
			xLblString = 'index [1/index_{max}]';
			
			switch prop
				case 'k'
% 					xx = obj.s;
					yy = obj.k;
% 					xLblString = 's [m]';
					yLblString = 'curvature [1/m]';
				
				case 'phi'
% 					xx = obj.s;
					yy = obj.phi;
% 					xLblString = 's [m]';
					yLblString = 'phi [rad]';
					
				case 's'
% 					n = length(obj.x);
% 					xx = (0:n-1)/n;
					yy = obj.s;
% 					xLblString = 'index [1/indMax]';
					yLblString = 's [m]';
					
				otherwise
					error(['Unknown property string. ',...
						'Type help plotdiff_ for valid property strings.']);
					
			end%switch
			
			% create SEGDAT object for plotting
			obj_new = segDat(...
				xx,...
				yy,...
				obj.s,...
				obj.k,...
				obj.phi,...
				obj.type,...
				obj.nbr);
			
			% plot
			h = plotdiff(obj_new,fh);
			
			% apply plot styles
			axis normal
			grid on
			xlabel(xLblString);
			ylabel(yLblString);
			
		end%fcn
		
		
		function h = plottangent(obj,ind,varargin) 
		%PLOTTANGENT	Plot the street segment and specified tangents.
		%	PLOTTANGENT(OBJ,IND) plots OBJ.y over OBJ.x, the tangents of
		%	the elements of indices IND and highlights the elements IND by
		%	a marker '*'.
		%	
		%	PLOTTANGENT(OBJ,IND,C) additionally applies the color C to the
		%	tangent and the marker.
		%	
		%	PLOTTANGENT(OBJ,IND,C1,...,CN) applies the colors C1,...,CN to
		%	IND(1),...,IND(N). If length(IND) > N, the pattern C1,...,CN
		%	will be extended to match length(IND). If length(IND) < N, the
		%	extra colors will be ignored.
		%	
		%	H = PLOTTANGENT(...) returns a handle H to lineseries objects.
		%	H(1,1) is the street segment-handle, H(i+1,1) and H(i+1,2) the
		%	marker- and the tangent-handle of IND(i) respectively.
		%	
		%	See also SEGDAT/PLOT, SEGDAT/PLOTDIFF.
		
			%%% handle input arguments
			% check class input ind
			if ~isnumeric(ind)
				error('Input argument IND has to be numeric.');
			end%if
			
			% check dimension of input ind
			if (~isvector(ind) && ~isempty(ind))
				error('Input argument IND has to be of size 0xN or Nx0.');
			end%if
			
			% apply plot options if unspecified
			if nargin < 3
				tangentColor = {'r'};
			else
				tangentColor = varargin;
			end%if
			
			% extend the color specification to the length of ind
			if length(ind) > length(tangentColor)
				tangentColor = repmat(tangentColor,1,ceil(length(ind)/length(tangentColor)));
			end
			
			% init handle to lineseries objects
			h = zeros(1+length(ind),2);
			
			% plot the segment and get corresponding axis limtis
			h(1,1) = plot(obj);
			h(1,2) = NaN;
			xLimits = xlim;
			yLimits = ylim;
			
			% check if some indexes are out of range
			if any(ind > length(obj.x))
				warning('MATLAB:lkaSegment:plottangent:indxOutOfRange',...
					['Some indexes exceed the number of elements ',...
					'of segment data which is ',num2str(length(obj.x)),'. ',...
					'Therefore their tangents will not be plotted.'])
				ind(ind > length(obj.x)) = [];
			end%if
			
			% plot the tangents and the corresponding segment-elements
			hold on
			for i = 1:length(ind)
				iind = ind(i);
				
				% marker of tangent point
				h(i+1,1) = plot(obj.x(iind),obj.y(iind),...
					'Marker','*','MarkerSize',7,'Color',tangentColor{i});
				
				% length of tangent
				[r1,r2] = obj.scaleTangentToAxis(xLimits,yLimits,...
					[obj.x(iind) obj.y(iind)],obj.phi(iind));
				
				% start/end point of tangent
				Pstart = [...
					obj.x(iind)+r2*cos(obj.phi(iind));...
					obj.y(iind)+r2*sin(obj.phi(iind))];
				Pstop = [...
					obj.x(iind)+r1*cos(obj.phi(iind));...
					obj.y(iind)+r1*sin(obj.phi(iind))];
				
				% draw the tangent using N arrows sticked together
				N = 25;
				xq = linspace(Pstart(1),Pstop(1),N);
				yq = linspace(Pstart(2),Pstop(2),N);
				% quiver would draw one arrow at every point of xq/yq of
				% length uq/vq; since the last element of xq/yq is the end
				% point, no arrow has to be drawn there
				xq = xq(1:end-1);
				yq = yq(1:end-1);
				uq = ones(size(xq))*(Pstop(1)-Pstart(1))/N;
				vq = ones(size(yq))*(Pstop(2)-Pstart(2))/N;
				scale = 0;
				h(i+1,2) = quiver(xq,yq,uq,vq,scale,'Color',tangentColor{i});
				
			end%for
			hold off
			
			% set the axis limtis corresponding to segment data
			axis([xLimits,yLimits]);
			
		end%fcn
		
		
		function plotLaneTracking(...
				sd_global,xyCG_global,psi_global,LAD,indx,...
				sd_T,xyCG_T)
			
			
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
			% desired path (global)
			h = plottangent(sd_global,indx);
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
				if ~isa(sd_T,'segDat')
					warning('SEGDAT:plotLaneTracking:class',...
						'Input argument SD_T not of class SEGDAT!');
% 					return;
					disp('Try to convert to class SEGDAT...')
					sd_T = segDat(...
						sd_T.x,...
						sd_T.y,...
						sd_T.s,...
						sd_T.k,...
						sd_T.phi,...
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
		
		function h = plot_raw(obj,varargin)
		%PLOT_RAW	Basic plot of street segment.
		%	To be used by public plot methods to avoid multiple calls to
		%	plot styles like AXIS, TITLE, XLABEL, ...!
		%
		% See also SEGDAT/PLOT.
		
			% get current status of 'NextPlot' property
			axh = gca;
			npState = get(axh,'NextPlot');
			
			% plot street segment
			h(1) = plot(obj.x,obj.y,varargin{:});
			
			% highlight first element as starting position
			set(axh,'NextPlot','add');
% 			h(2) = plot(obj.x(1),obj.y(1),varargin{:},'Marker','o');
% 			plot(obj.x(1),obj.y(1),varargin{:},'Marker','o');
			color = get(h(1),'Color');
			plot(obj.x(1),obj.y(1),varargin{:},...
				'Color',color,...
				'MarkerFaceColor',color);
			set(axh,'NextPlot',npState); % reset to initial state
			
		end%fcn
		
		
		function h = quiver_raw(obj,varargin)
		%QUIVER_RAW		Basic quiver plot of street segment.
		%	To be used by public plot methods to avoid multiple calls to
		%	plot styles like AXIS, TITLE, XLABEL, ...!
		%	
		% See also SEGDAT/QUIVER.
		
			% plot street segment
			u = diff(obj.x);
			v = diff(obj.y);
			h(1) = quiver(obj.x(1:end-1),obj.y(1:end-1),u,v,...
				varargin{:},'AutoScale','off');
			
		end%fcn
		
		
		function cellStr = getLegendCellString(obj)
			
			if all(obj.type(1) == obj.type)
				gettype = obj.curvTypes{obj.type(1) + 2};
			else
				gettype = 'mixed';
			end%if
			lengthStr	= [sprintf('%.2f',obj.s(end)),' m'];
			pointsStr	= sprintf('%.0d',length(obj.s));
			startStr	= sprintf('(%g;%g)',obj.x(1),obj.y(1));
			endStr		= sprintf('(%g;%g)',obj.x(end),obj.y(end));
			cellStr		= {...
				['Street segment of type ''',gettype,''''],...
				['length: ',lengthStr,' using ',pointsStr,' points'],...
				['moves from ',startStr,' \rightarrow ', endStr]};
			
		end%fcn
		
	end%methods
	
	
	%%% SET-Methods
	methods
		
		function obj = set.type(obj,value)
			
			if any(value < -1) || any(value > 3)
				error(['Street segment TYPE out of range! ',...
					'Type help SEGDAT/TYPE for valid values.']);
			else
				obj.type = value;
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
			
			%{ 
			this approach fails since some angles alphX do not exist if the
			tangents point is located in one of the plot corners.
			
			% map phi to [0,2*pi)
			phi = mod(phi,2*pi);
			
			% get the angles from tangent point to XLIMITS/YLIMITS corners
			alph1 = mod(atan2(yMax-yT,xMax-xT)+2*pi,2*pi);
			alph2 = mod(atan2(yMax-yT,xMin-xT)+2*pi,2*pi);
			alph3 = mod(atan2(yMin-yT,xMin-xT)+2*pi,2*pi);
			alph4 = mod(atan2(yMin-yT,xMax-xT)+2*pi,2*pi);
			
			% calc r1
			if 0 <= phi && phi <= alph1
				r1 = (xMax-xT)/cos(phi);
				
			elseif alph1 < phi && phi <= alph2
				r1 = (yMax-yT)/sin(phi);
				
			elseif alph2 < phi && phi <= alph3
				r1 = (xMin-xT)/cos(phi);
				
			elseif alph3 < phi && phi <= alph4
				r1 = (yMin-yT)/sin(phi);
				
			elseif alph4 < phi && phi <= 2*pi
				r1 = (xMax-xT)/cos(phi);
				
			else
				error('Fatal error at calculating r1!');
			end%if
			
			% calc r2
			phi_ = mod(phi + pi,2*pi);
			if 0 <= phi_ && phi_ <= alph1
				r2 = (xMax-xT)/cos(phi);
				
			elseif alph1 < phi_ && phi_ <= alph2
				r2 = (yMax-yT)/sin(phi);
				
			elseif alph2 < phi_ && phi_ <= alph3
				r2 = (xMin-xT)/cos(phi);
				
			elseif alph3 < phi_ && phi_ <= alph4
				r2 = (yMin-yT)/sin(phi);
				
			elseif alph4 < phi_ && phi_ <= 2*pi
				r2 = (xMax-xT)/cos(phi);
				
			else
				error('Fatal error at calculating r2!');
			end%if
			%}
			
			
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
		
		
		function fid = openFileWithExt(fn,fExt_set)
			
			if nargin > 1
				[~,fName,fExt] = fileparts(fn);
				if isempty(fExt)
					warning('segDat:write2file:fileExtension',...
						'Adding file extension ''%s'' to file name!',...
						fExt_set);
				elseif ~strcmp(fExt,fExt_set)
					warning('segDat:write2file:fileExtension',...
						'Replacing file extension ''%s'' by ''%s''!',...
						fExt,fExt_set);
				end%if
				fn = [fName,fExt_set];
			end%if
			
			% open file with write-permission
			fid = fopen(fn,'w');
				
		end%fcn
		
		
		function sd = xy2segDat(x,y)
		%XY2SEGDAT	Convert x/y coordinates to SEGDAT class. 
		%   OBJ = XY2SEGDAT(X,Y) creates the instace OBJ of class SEGDAT
		%   from X/Y-coordinates given in meter. X and Y must be vectors
		%   with the same number of elements!
		% 
		%   OBJ = XY2SEGDAT(XY) X/Y-coordinates are given in the array XY
		%   of size 2-by-n.
		%
		%	Property S is calculated using cumulative sum.
		%	Properties K and PHI are set to NaN.
		%	Property TYPE is set to -1.
		%	Property NBR is set to 1.
		
		
			narginchk(1,2);
		
			%%% handle input arguments
			if nargin < 2
				if size(x,1) ~= 2
					error('When using one argument syntax, array must have 2 rows!');
				else
					y = x(2,:);
					x = x(1,:);
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
			
			
			%%% convert to SEGDAT
			sd = segDat(...
				x,...
				y,...	
				[0,cumsum(sqrt(diff(x).^2 + diff(y).^2))],...
				NaN(size(x)),...
				NaN(size(x)),...
				-1,...
				1);
			
		end%fcn

	end%methods
	
	
	%%% Test-Methods
	methods (Static, Hidden)
		
		function test_changeSignOfCurvature(obj)
			
			if nargin < 1
				b = lkaSegmentCircle([],5/3*pi,4/2*pi,50);
				b = shift(b,[20,10]);
				obj = b.segmentData;
			end%if
			
			fig = figure;
			
			ind = [1,length(obj.x)];
			h = plottangent(obj,ind);
			set(h(1,1),'Color','b','LineWidth',1);
			set(h(2:end,1),'Color','c','Marker','o','MarkerFaceColor','c');
			set(h(2:end,2),'Color','c');
			
			obj_ = changeSignOfCurvature(obj);
			hold on
			h_ = plottangent(obj_,ind);
			set(h_(1,1),'Color','r','LineWidth',1);
			set(h_(2:end,1),'Color','m','Marker','h','MarkerFaceColor','m');
			set(h_(2:end,2),'Color','m');
			
			axis auto
			legend([h(1,1),h_(1,1)],'original','curvature*(-1)','location','Best')
			pause
			close(fig)
			
		end%fcn
		
		
		function test_rotate(P0,phi)
			
			if nargin < 2; phi = pi/2; end%if
			if nargin < 1; P0 = [20 0]; end%if
			
			fig = figure;
			b = lkaSegmentCircle([],3/2*pi,4/2*pi,50);
			
			sd_b = b.segmentData;
			
			sd_b = shiftTo(sd_b,P0);
			
			plotdiff(sd_b);
			
			plotdiff(rotate(sd_b,phi));
			
			pause
			close(fig)
			
		end%fcn
		
		
		function test_reverseDirection(obj)
			
			fig = figure;
			
			obj_ = reverseDirection(obj);
			
			ind = 1:10:21;
			h = plottangent(obj,ind,'b');
			
			hold on
			h = plottangent(obj_,ind,'r');
			set(h(1,1),'LineStyle',':','Color','r');
			
			pause
			close(fig)
			
		end%fcn
		
		
		function indFailed = test_scaleTangentToAxis(nbr)
			
			
			% define some test cases and their desired results
			tc = [...
				%solution	xLimits		yLimits		point		angle
				%
				% tangent point at lower left corner
				{[+10 0],	[0 10],		[0 10],		[0 0],		0};...
				{[+10 0],	[0 10],		[0 10],		[0 0],		pi/2};...
				{[0 -10],	[0 10],		[0 10],		[0 0],		pi};...
				{[0 -10],	[0 10],		[0 10],		[0 0],		3*pi/2};...
				%
				% tangent point at lower right corner
				{[0 -10],	[0 10],		[0 10],		[10 0],		0};...
				{[+10 0],	[0 10],		[0 10],		[10 0],		pi/2};...
				{[+10 0],	[0 10],		[0 10],		[10 0],		pi};...
				{[0 -10],	[0 10],		[0 10],		[10 0],		3*pi/2};...
				%
				% tangent point at upper right corner
				{[0 -10],	[0 10],		[0 10],		[10 10],	0};...
				{[0 -10],	[0 10],		[0 10],		[10 10],	pi/2};...
				{[+10 0],	[0 10],		[0 10],		[10 10],	pi};...
				{[+10 0],	[0 10],		[0 10],		[10 10],	3*pi/2};...
				%
				% tangent point at upper left corner
				{[+10 0],	[0 10],		[0 10],		[0 10],		0};...
				{[0 -10],	[0 10],		[0 10],		[0 10],		pi/2};...
				{[0 -10],	[0 10],		[0 10],		[0 10],		pi};...
				{[+10 0],	[0 10],		[0 10],		[0 10],		3*pi/2};...
				%
				% tangent point at upper left corner
				% k = tan(phi)
				% d = k*x - y(x) = 5 - 5*tan(phi)
				% r1 = sqrt(5^2 + (5-d)^2)
				% r2 = -r1
				{[sqrt(25+25*tan(pi/8)^2) -sqrt(25+25*tan(pi/8)^2)],...
							[0 10],		[0 10],		[5 5],		pi/8};
			];
		
% 			cell2struct(tc,{'sol','xLimits','yLimits','point','angle'},2);
			
			% handle input arguments
			if nargin < 1
				nbr = 1:size(tc,1);
			end%if
			
			% select test cases to perform
			tc = tc(nbr,:);
			nbrOfTestCases = size(tc,1);
			
			% pre-allocation
			isTestSuccesfull = false(nbrOfTestCases,1);
			
			% run test cases
			for i = 1:nbrOfTestCases
				
				[r1,r2] = segDat.scaleTangentToAxis(...
					tc{i,2},...
					tc{i,3},...
					tc{i,4},...
					tc{i,5});
				
				isTestSuccesfull(i) = all(tc{i,1} == [r1 r2]);
				
			end%for
			
			% evaluation and command line output
			nbrSuccesfull = sum(isTestSuccesfull);
			nbrFailed = sum(~isTestSuccesfull);
			fprintf('   Test results: \n')
			fprintf('     succesfull: %3d of %d.\n',nbrSuccesfull,nbrOfTestCases)
			fprintf('     failed:     %3d of %d.\n',nbrFailed,nbrOfTestCases)
			
			ind = 1:length(tc);
			indFailed = ind(~isTestSuccesfull);
			formatString = repmat('%d, ',1,nbrFailed);
			formatString = [formatString(1:end-2),'.'];
			if ~isempty(indFailed)
				fprintf(['     test cases failed: ',formatString,'\n'],...
					indFailed);
			end%if
			fprintf('\n');
			
		end%fcn
		
		
		function test_shiftTo()
			
			fig = figure;
			a = lkaSegmentStraight([],100,pi/8);
			b = lkaSegmentCircle([],3/2*pi,4/2*pi,50);
			c = lkaSegmentClothoid([],1/100,0,0,100);
			
			sd_a = a.segmentData;
			sd_b = b.segmentData;
			sd_c = c.segmentData;
			
			sd_a = shiftTo(sd_a,[10 30]);
			sd_b = shiftTo(sd_b,[20 20]);
			sd_c = shiftTo(sd_c,[30 10]);
			
			plotdiff(sd_a);
			plotdiff(sd_b);
			plotdiff(sd_c);
			
			plotdiff(shiftTo(sd_a));
			plotdiff(shiftTo(sd_b));
			plotdiff(shiftTo(sd_c));
			
			pause
			close(fig)
			
		end%fcn
		
		
		function test_plotdiff(sd)
			
			if nargin < 1
				a = lkaSegmentStraight([],200,0);
				b = lkaSegmentClothoid([],0,0.01,a.segmentData.phi(end),200);
				c = lkaSegmentCircle([],b.segmentData.phi(end)-pi/2,pi*1/2,300);
				d = lkaSegmentClothoid([],0.01,0,c.segmentData.phi(end),200);
				e = lkaSegmentStraight([],200,d.segmentData.phi(end));
				sd = a + b + c + d + e;
			end%if
			
			plotdiff(sd);
			
		end%fcn
		
		
		function test_plotdiff_(sd)
			
			if nargin < 1
				a = lkaSegmentStraight([],200,0);
				b = lkaSegmentCircle([],-pi/2,pi/2,500);
				c = lkaSegmentStraight([],200,pi);
				d = lkaSegmentCircle([],pi/2,3*pi/2,500);
				sd = a + b + c + d;
			end%if
			
			% x/y
			subplot(4,1,1)
			plotdiff(sd);
			
			% curve length
			subplot(4,1,2)
			plotdiff_(sd.segmentData,[],'s');
			title('')
			
			% curvature
			subplot(4,1,3)
			plotdiff_(sd.segmentData,[],'k');
			title('')
			
			% tangent angle
			subplot(4,1,4)
			plotdiff_(sd.segmentData,[],'phi');
			title('')
			
		end%fcn
		
	end%methods
	
end%class
