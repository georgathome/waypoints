function cellStr = getTitleCellString(obj)

	if numel(obj) > 1
		cellStr = {sprintf('Comparing %d paths', numel(obj))};
	else

		if all(obj.Type(1) == obj.Type)
			typeStr = obj.CurvTypes{obj.Type(1) + 2};
			typeStr(1) = upper(typeStr(1));
		else
			typeStr = 'Undefined Waypoints';
		end%if
		cellStr	= {...
			sprintf('%s (N=%u;  s=%.2f)', ...
			typeStr, numwp(obj), obj.s(end)),...
			sprintf(' (%.2g;%.2g) \\rightarrow (%.2g;%.2g)', ...
			obj.x(1), obj.y(1), obj.x(end), obj.y(end)),...
			};

	end%if

end%fcn
