function test_changeSignOfCurvature(obj)
			
	if nargin < 1
		b = LkPathCircle([],5/3*pi,4/2*pi,50);
		b = shiftBy(b,[20,10]);
		obj = b.pathData;
	end%if

	figure;
 
	ind = [1,length(obj.x)];
	h = plottangent(obj,ind);
	set(h(1,1),'Color','b','LineWidth',1);
	set(h(2:end,1),'Color','c','Marker','o','MarkerFaceColor','c');
	set(h(2:end,2),'Color','c');

	obj_ = changeSignOfCurvature(obj);
	hold on
	h_ = plottangent(obj_,ind);
	set(h_(1,1),'Color','r','LineWidth',1);
	set(h_(2:end,1),'Color','m','Marker','h','MarkerFaceColor','m');
	set(h_(2:end,2),'Color','m');

	axis auto
	legend([h(1,1),h_(1,1)],'original','curvature*(-1)','location','Best')

end%fcn
