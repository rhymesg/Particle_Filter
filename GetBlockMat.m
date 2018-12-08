function [ ele ] = GetBlockMat( Mat_in, i, j, r, c )
%GETMAT 이 함수의 요약 설명 위치
%   자세한 설명 위치

ele = Mat_in(r*(i-1)+1:r*i, c*(j-1)+1:c*j);

end

