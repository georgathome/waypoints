function [head, curv] = curvHeadFrom_x_y(x, y)

g1_x = gradient(x);
g1_y = gradient(y);
g2_x = gradient(g1_x);
g2_y = gradient(g1_y);

head = unwrap(atan2(g1_y, g1_x));
if nargout > 1
	curv = (g1_x.*g2_y - g2_x.*g1_y) ./ (g1_x.^2 + g1_y.^2).^1.5;
end

end%fcn
