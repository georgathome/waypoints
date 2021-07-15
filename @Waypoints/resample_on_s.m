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
