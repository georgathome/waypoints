function s = sFrom_x_y(x, y)

s = [0, cumsum(sqrt(diff(x).^2 + diff(y).^2))];

end%fcn
