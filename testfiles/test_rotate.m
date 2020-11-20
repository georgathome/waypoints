function test_rotate(P0, phi)
			
	if nargin < 2; phi = pi/2; end%if
	if nargin < 1; P0 = [20 0]; end%if

	fig = figure;
	b = LkPathCircle([],3/2*pi,4/2*pi,50);

	sd_b = b.pathData;

	sd_b = shiftTo(sd_b,P0);

	plotdiff(sd_b);

	plotdiff(rotate(sd_b,phi));

	pause
	close(fig)

end%fcn
