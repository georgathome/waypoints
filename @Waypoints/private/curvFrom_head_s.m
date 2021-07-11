function curv = curvBy_head_s(head, s)

curv = gradient(head)./gradient(s);

end%fcn
