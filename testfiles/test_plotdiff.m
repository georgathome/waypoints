function h = test_plotdiff(obj)
			
	if nargin < 1
		a = LkPathStraight([], 200, 0);
		b = LkPathClothoid([], 0, 0.01, a.pathData.head(end), 200);
		c = LkPathCircle([], b.pathData.head(end)-pi/2, pi*1/2, 300);
		d = LkPathClothoid([], 0.01, 0, c.pathData.head(end), 200);
		e = LkPathStraight([], 200, d.pathData.head(end));
		obj = a + b + c + d + e;
	end%if

	h = plotdiff(obj);

end%fcn