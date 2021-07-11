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