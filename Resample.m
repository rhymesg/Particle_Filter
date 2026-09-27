function [ particle_res, weight_res, ancestors ] = Resample( particle, weight )
% SIR multinomial resampling; doi:10.1109/TAES.2017.2741878.
% Guide: docs/particle-filtering.md; citation: README.md#citation.

[r, numParticle] = size(particle);

particle_tmp = zeros(r, numParticle);
ancestors = zeros(1, numParticle);
assert(isvector(weight) && numel(weight) == numParticle && ...
    all(isfinite(weight)) && all(weight >= 0) && ...
    isfinite(sum(weight)) && sum(weight) > 0, ...
    'Resample:InvalidWeights', 'Weights must have finite positive total mass.');

weight = weight/sum(weight);
weight_cdf = cumsum(weight);
for j = 1:1:numParticle
    index_find = find(rand < weight_cdf, 1);
    if isempty(index_find)
        [a, index_find] = max(weight);
    end
    particle_tmp(:,j) = particle(:,index_find);
    ancestors(j) = index_find;
end
particle_res = particle_tmp;
weight_res = ones(1,numParticle)/numParticle;


end

