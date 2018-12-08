function [ particle_res, oosmSucceed ] = OOSM( particle_mi, weight_un, dcov_pl, skipped, sig_meas, DEM, oosmSucceed )
%OOSM 이 함수의 요약 설명 위치
%   자세한 설명 위치

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

