function [ Mat_out ] = SetBlockMat( Mat_in, ele, i, j, r, c )
%SetMat 이 함수의 요약 설명 위치
%   자세한 설명 위치

Mat_out = Mat_in;

Mat_out(r*(i-1)+1:r*i, c*(j-1)+1:c*j) = ele;

end

