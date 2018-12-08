function [ out ] = covAnal( cov )
%COVANAL 이 함수의 요약 설명 위치
%   자세한 설명 위치

s_cov = sqrt(cov);

% Det = det(s_cov);

dstd = sqrt(s_cov(1,1)^2 + s_cov(2,2)^2);

stdtr = sqrt(trace(cov));

out = [dstd; stdtr];


end

