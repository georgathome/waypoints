function cellStr = getTitleCellString(obj)

	if numel(obj) > 1
		cellStr = {sprintf('Comparing %d Waypoints', numel(obj))};
	else
		cellStr	= {...
			sprintf('Waypoints: N=%u; s=%.1f; (%.2g;%.2g) \\rightarrow (%.2g;%.2g)', ...
			numwp(obj), obj.s(end), obj.x(1), obj.y(1), obj.x(end), obj.y(end)),...
			};
	end%if

end%fcn
