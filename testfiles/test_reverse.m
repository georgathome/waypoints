function test_reverse(obj)
			
	fig = figure;

	ind = 1:10:21;
	obj.Name = 'Forward';
	h = plottangent(obj, ind, 'Color','b');
	set(h(1,1), 'Marker','none');

	hold on
	obj = reverse(obj);
	obj.Name = 'Reverse';
	h = plottangent(obj, ind, 'Color','r');
	set(h(1,1), 'Marker','none');
	hold off

	legend('show');
	pause
	close(fig)

end%fcn