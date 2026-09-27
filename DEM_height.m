function [ height ] = DEM_height( pos, DEM )
% Terrain-height interpolation for TRN; doi:10.1109/TAES.2017.2741878.
% Guide: docs/simulation.md; citation: README.md#citation.

resolution = DEM.resolution;

[r, c] = size(DEM.DB);

indx = floor(pos(1)/resolution);
indy = floor(pos(2)/resolution);

X = indx-2:1:indx+2;
Y = indy-2:1:indy+2;

DB_part = DEM.DB(indx-2:indx+2, indy-2:indy+2);

height = interp2(Y, X, DB_part, pos(2)/resolution, pos(1)/resolution, 'cubic'); % x, y ¹Ý´ë·Î

end

