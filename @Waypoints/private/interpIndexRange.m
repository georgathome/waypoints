function [indl,indu,ind] = interpIndexRange(ind,indMinMax,m)
%INTERPINDEXRANGE	Index range for interpolation.
%	[INDL,INDU] = INTERPINDEXRANGE(IND,INDMINMAX,M) returns the lower and
%	upper indices INDL and INDU within index boundaries [INDMINMAX(1)
%	INDMINMAX(2)] using an index difference M for given indices IND.
%	
%	Indices INDL/INDU define the range of interpolation, basically IND-M
%	and IND+M but ensure INDL > INDMINMAX(1) and INDU < INDMINMAX(2).


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
