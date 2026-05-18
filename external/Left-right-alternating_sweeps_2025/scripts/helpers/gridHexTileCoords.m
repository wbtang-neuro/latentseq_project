function [x, y] = gridHexTileCoords(gridAxes)

f = tan(pi/6);
tileVertices = rotate2d(gridAxes*f, deg2rad(30));

x = [];
y = [];
c = 0;

for pol = [-1, 1]
    for a = 1:3
        c = c+1;
        x(c) = tileVertices(a, 1)*pol;
        y(c) = tileVertices(a, 2)*pol;
    end
end

end