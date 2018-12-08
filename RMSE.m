function [ X_err_RMS, Y_err_RMS, d_err_RMS ] = RMSE( x_err, length, numMonte )
%RMSE 이 함수의 요약 설명 위치
%   자세한 설명 위치

X_err = zeros(length, numMonte);
Y_err = zeros(length, numMonte);
d_err = zeros(length, numMonte);

for k = 1:1:length
    for i = 1:1:numMonte
        d_err(k,i) = norm(x_err(:,k,i)); 
    end
end
X_err(:,:) = x_err(1,:,:);
Y_err(:,:) = x_err(2,:,:);

X_err_RMS = sqrt(mean(X_err.^2,2));
Y_err_RMS = sqrt(mean(Y_err.^2,2));
d_err_RMS = sqrt(mean(d_err.^2,2));

end

