function write2file_CarMaker(fn, obj, varargin)

	% open file
	[fid, fileExt] = openFileWithExt(fn);

	switch fileExt
		case '.csv'
			% x .. x-coordinate [m]
			% y .. y-coordinate [m]
			% z .. altitude [m]
			% q .. slope
			% wl/wr .. track width left/right [m]
			% ml/mr .. margin width left/right [m]
			N = numwp(obj);
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
				[obj.x, obj.y, zeros(N,1), zeros(N,1), ...
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
