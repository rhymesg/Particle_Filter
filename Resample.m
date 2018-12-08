function [ particle_res, weight_res ] = Resample( particle, weight )
%RESAMPLE 이 함수의 요약 설명 위치
%   자세한 설명 위치

[r, numParticle] = size(particle);

particle_tmp = zeros(r, numParticle);

weight = weight/sum(weight);
weight_cdf = cumsum(weight);
for j = 1:1:numParticle;
    index_find = find(rand <= weight_cdf, 1);
    if isempty(index_find)
        [a, index_find] = max(weight);
    end
    particle_tmp(:,j) = particle(:,index_find);
end
particle_res = particle_tmp;
weight_res = ones(1,numParticle)/numParticle;


end

