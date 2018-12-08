function [ cluster ] = Cluster( particle, center )
%CLUSTER 이 함수의 요약 설명 위치
%   자세한 설명 위치

numParticle = size(particle,2);
numMode = size(center,2);
dists = zeros(1,numMode);
cluster = cell(numMode, 1);
for n = 1:1:numParticle
    for m = 1:1:numMode
        dists(m) = norm (particle(:,n) - center(:,m));
    end
    
    [val, ind] = min(dists);
    
    cluster{ind} = [cluster{ind}, particle(:,n)];
    
end




end

