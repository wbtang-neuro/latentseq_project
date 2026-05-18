function pos = posDistShift(pos, direction, distance)
dx = cos(direction) .* distance;
dy = sin(direction) .* distance;
pos = pos + [dx, dy];
end