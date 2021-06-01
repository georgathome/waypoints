function s = sFrom_x_y(x, y)

s = [0; cumsum(hypot(diff(x), diff(y)))];

end%fcn
