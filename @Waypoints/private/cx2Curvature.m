function curv = cx2Curvature(d1x, d1y, d2x, d2y)
curv = (d1x.*d2y - d2x.*d1y) ./ (d1x.^2 + d1y.^2).^1.5;
end%fcn
