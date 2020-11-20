function p = test_resample()

	% Create analytical test data
	a = LkPathStraight(5, 100, 0);
	b = LkPathCircle(2, -pi/2, 0, 50);
	c = LkPathClothoid(20, 0, 1/15, 90*pi/180, 70);

	p1 = a.pathData + b.pathData + c.pathData;

	ds = 3;
	p2.Name = 'linear';
	p2 = resample(p1, ds, p2.Name);
	p3.Name = 'spline';
	p3 = resample(p1, ds, p3.Name);


	ax = gobjects(4, 1);
	ax(1) = subplot(3,2,[1,3,5]);
	h1 = plotdiff(p1);
	hold on
	h2 = plotdiff(p2);
	h3 = plotdiff(p3);
	hold off
	set(h1, 'Marker','o');
	set(h2, 'Marker','d');
	set(h3, 'Marker','.')
	title('Test: resample()');
	legend(...
		'original straight','original circle','original clothoid',...
		'linear straight','linear circle','linear clothoid',...
		'spline straight','spline circle','spline clothoid');

	ax(2) = subplot(3,2,2);
	test_resample_sub(p1, p2, p3, 'head');
	legend('original','linear','spline');

	ax(3) = subplot(3,2,4);
	test_resample_sub(p1, p2, p3, 'curv');

	ax(4) = subplot(3,2,6);
	test_resample_sub(p1, p2, p3, 'nbr');

	linkaxes(ax(2:4),'x');

	p = [p1; p2; p3];
	
end%fcn

function test_resample_sub(obj0, obj1, obj2, field)
	plot(obj0.s, obj0.(field), 'ro');
	hold on
	plot(obj1.s, obj1.(field), 'b.');
	plot(obj2.s, obj2.(field), 'g.');
	hold off
	grid on
	xlabel('s [m]');
	title(field);
end%fcn
