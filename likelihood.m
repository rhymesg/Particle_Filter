function [ p ] = likelihood( Z_est, Z_mea, sig_meas )
%LIKELIHOOD 이 함수의 요약 설명 위치
%   자세한 설명 위치

    var = sig_meas^2;
    p = (1/sqrt(2*pi*var)) * exp( -(Z_mea - Z_est)^2/(2*var) );

end

