function [ height ] = DEM_height( pos, DEM )
%DTED_HEIGHT 이 함수의 요약 설명 위치
%   자세한 설명 위치

resolution = DEM.resolution;

[r, c] = size(DEM.DB);

indx = floor(pos(1)/resolution);
indy = floor(pos(2)/resolution);

X = indx-2:1:indx+2;
Y = indy-2:1:indy+2;

DB_part = DEM.DB(indx-2:indx+2, indy-2:indy+2);

height = interp2(Y, X, DB_part, pos(2)/resolution, pos(1)/resolution, 'cubic'); % x, y 반대로

end

