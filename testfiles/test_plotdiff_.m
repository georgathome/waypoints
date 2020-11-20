function h = test_plotdiff_(obj)
			
	if nargin < 1
		a = LkPathStraight([], 200, 0);
		b = LkPathCircle([], -pi/2, pi/2, 500);
		c = LkPathStraight([], 200, pi);
		d = LkPathCircle([], pi/2, 3*pi/2, 500);
		obj = a.pathData + b.pathData + c.pathData + d.pathData;
	end%if

	ax = gobjects(4, 1);
	h = gobjects(4, max(obj.nbr));
	% x/y
	ax(1) = subplot(3,2,[1,3,5]);
	h(1,:) = plotdiff(obj);

	% curve length
	ax(2) = subplot(3,2,2);
	h(2,:) = plotdiff_(obj, [], 's');
	title('');
	xlabel('');

	% tangent angle
	ax(3) = subplot(3,2,4);
	h(3,:) = plotdiff_(obj, [], 'head');
	title('');
	xlabel('');

	% curvature
	ax(4) = subplot(3,2,6);
	h(4,:) = plotdiff_(obj, [], 'curv');
	title('');

	linkaxes(ax(2:end), 'x');

end%fcn