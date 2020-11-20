function test_shiftTo()
			
	a = LkPathStraight([], 100, pi/8);
	b = LkPathCircle([], 3/2*pi, 4/2*pi, 50);
	c = LkPathClothoid([], 1/100, 0, 0, 100);

	sd_a = a.pathData;
	sd_b = b.pathData;
	sd_c = c.pathData;

	sd_a = shiftTo(sd_a, [10 30]);
	sd_b = shiftTo(sd_b, [20 20]);
	sd_c = shiftTo(sd_c, [30 10]);

	fig = figure;
	hold on
	plotdiff(sd_a);
	plotdiff(sd_b);
	plotdiff(sd_c);

	plotdiff(shiftTo(sd_a));
	plotdiff(shiftTo(sd_b));
	plotdiff(shiftTo(sd_c));
	hold off

	pause
	close(fig)

end%fcn