function [ particle_res, oosmSucceed ] = OOSM( particle_mi, weight_un, dcov_pl, skipped, sig_meas, DEM, oosmSucceed )
% Stored-particle OOSM variant; doi:10.1109/TAES.2017.2741878.
% Guide: docs/particle-filtering.md; citation: README.md#citation.

[a M] = size(skipped);
[r, numParticle] = size(particle_mi);

for m = 1:1:M;
    oosmindex = skipped(m).index;
    
    
    weight_tmp = weight_un;
    particle_oosm = skipped(m).particle;
    for n = 1:1:numParticle
        z_est = DEM_height(particle_oosm(:,n),DEM);
        weight_un(n) = weight_un(n) * likelihood(z_est, skipped(m).z, sig_meas);
    end

    weight = weight_un/sum(weight_un);
    [particle_star, w] = Resample(particle_mi, weight);
    
    dcov_oosm = det(cov(particle_star'));
    
end

particle_res = particle_star;


end

