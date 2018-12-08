%%
% Author: Youngjoo Kim
% Please cite the following paper if you find this code helpful:
% Youngjoo Kim et al., "Utilizing Out-of-Sequnece Measurement for
%   Ambiguous Update in Particle Filtering", IEEE Transactions on Aerospace
%   and Electronic Systems, 54(1), 2018.

%% setting
load('DB_part.mat')

% switches (1: run, 0: don't run)
RUN_DATA = 1;   % data loading
RUN_OOSM = 1;   % proposed method with out-of-sequence measurement update
RUN_PF = 1;     % standard particle filter
RUN_MPF = 1;    % mixture particle filter
RUN_APF = 1;    % auxiliary particle filter

simTime = 150;
dt = 1;
numParticle = 500;
numMonte = 1; % number of Monte-Carlo runs
length = simTime/dt + 1; % t = 0 at k = 1
time = 0:dt:simTime;

x_init = [2600; 2800];
x_init_err = [30; 30];

% sensor setting
sig_z = 15 + 4.71;  % std of actual measurement noise
sig_v = 2;          % std of actual process noise
sig_bias_z = 0;     % actual measurement noise bias
sig_bias_v = 0;     % actual process noise bias

% filter setting
sig_init = 40;      % std for initial P
sig_proc = 8;       % std for Q
sig_meas = 25;      % std for R

%% true trajectory & sensor data

if (RUN_DATA == 1)

spd = 42;
tlength = 300 + 1; % trajectory length
t_circle = round(560*pi/spd)+1; % time to make half circle

x_true = zeros(2,tlength);
v_true = zeros(2,tlength-1);
x_true(:,1) = x_init;
v_true(:,1) = [spd; 0];
z_meas = zeros(length, numMonte);
v_meas = zeros(2, length, numMonte);

for k = 2:1:tlength
    if (k <= 130)
        v_true(:,k-1) = [spd; 0];
    elseif (k <= 130 + t_circle)
        v_true(:,k-1) = spd*[cos(pi*(k-130)/t_circle); sin(pi*(k-130)/t_circle)];
    else
        v_true(:,k-1) = [-spd; 0];
    end
    
    x_true(:,k) = x_true(:,k-1) + v_true(:,k-1)*dt;
end

% figure;
% plot(x_true(1,:), x_true(2,:));
% axis equal;

bias_z = sig_bias_z*randn(1,numMonte);
bias_v = sig_bias_v*randn(2,numMonte);
for i = 1:1:numMonte
    rng('shuffle');
    for k = 2:1:length
        z_meas(k,i) = DEM_height(x_true(:,k), DEM) + sig_z*randn + bias_z(i);
        v_meas(:,k-1,i) = v_true(:,k-1) + sig_v*randn(2,1) + bias_v(:,i);
    end
end

end

%% simulation

%% OOSM
if (RUN_OOSM == 1)

x_est = zeros(2,length,numMonte);
x_err = zeros(2,length,numMonte);
dcov_res = zeros(2,length,numMonte);
particle = zeros(2,numParticle);
particle_mi = zeros(2,numParticle);
particle_pl = zeros(2,numParticle);
covIncrease = zeros(1,length);
oosmSucceed = zeros(1,length);
numMode = zeros(1,length);
disp(['-- OOSM started --']);
for i = 1:1:numMonte
    rng('shuffle');
    
    x_est(:,1,i) = x_init + x_init_err;
    x_err(:,1,i) = x_init_err;
    particle = [x_est(1,1,i)*ones(1,numParticle); x_est(2,1,i)*ones(1,numParticle)] + sig_init*randn(2,numParticle);
    weight = ones(1,numParticle)/numParticle;
    skipped = [];
    last.particle = particle;
    last.index = 1;
    for k = 2:1:length
        
        %%% time update
        for n = 1:1:numParticle
            particle(:,n) = particle(:,n) + v_meas(:,k-1,i)*dt;
            particle(:,n) = particle(:,n) + sig_proc*randn(2,1)*dt;
        end
        
        particle_mi = particle;
        dcov_mi = det(cov(particle_mi'));
        
        %%% measurement update
        for n = 1:1:numParticle       
            z_est = DEM_height(particle(:,n),DEM);
            weight(n) = weight(n) * likelihood(z_est, z_meas(k,i), sig_meas);
        end
        weight_un = weight;
        weight = weight/sum(weight);
        
        %%% resampling (SIR)
        [particle, weight] = Resample(particle, weight);
          
        particle_pl = particle;
        cov_pl = cov(particle_pl');
        dcov_pl = det(cov_pl);

%         numMode(k) = NumMode(particle_pl);
%         numMode(k);
        
        %%% cov test
        if (dcov_pl > dcov_mi) % if posterior cov > prior cov
            covIncrease(k) = covIncrease(k) + 1;
            
            particle = particle_mi; % skip update
            
            oosm.index = k;
            oosm.z = z_meas(k,i);
            oosm.v = v_meas(:,k-1,i);
            oosm.particle = particle_mi;
            skipped = [skipped oosm];
        else
            if (size(skipped,2) > 0) % if some measurements have been skipped
                [particle, oosmSucceed] = OOSM(particle_mi, weight_un, dcov_pl, skipped, sig_meas, DEM, oosmSucceed);
                skipped = [];
            end
            
            last.particle = particle; % store last successful update
            last.index = k;
        end
  
        dcov_res(:,k,i) = covAnal(cov(particle'));
        
        x_est(:,k,i) = mean(particle'); % MLSE estimate
        x_err(:,k,i) = x_est(:,k,i) - x_true(:,k);

    end
%     if (mod(i,1) == 0)
        disp([' OOSM Monte-Carlo step ' int2str(i) '/' int2str(numMonte) ' completed']);
%     end
end
disp(['-- OOSM completed --']);
x_est_OOSM = x_est;
x_err_OOSM = x_err;
particle_OOSM = particle;
covIncrease_OOSM = covIncrease;
numMode = numMode/numMonte;
dcov_res_OOSM = dcov_res;

[X_err_RMS_OOSM, Y_err_RMS_OOSM, d_err_RMS_OOSM] = RMSE(x_err_OOSM, length, numMonte);

end

%% PF
if (RUN_PF == 1)
x_est = zeros(2,length,numMonte);
x_err = zeros(2,length,numMonte);
particle = zeros(2,numParticle);
dcov_res = zeros(2,length,numMonte);
% dcov_mi = zeros(length,numMonte);
% dcov_pl = zeros(length,numMonte);
numMode = zeros(1,length);
signi = zeros(1,length);
gradx = zeros(1,length);
grady = zeros(1,length);
covIncrease = zeros(1,length);
dist = zeros(1,length);
disp(['-- PF started --']);
for i = 1:1:numMonte
    rng('shuffle');
    
    x_est(:,1,i) = x_init + x_init_err;
    x_err(:,1,i) = x_init_err;
    particle = [x_est(1,1,i)*ones(1,numParticle); x_est(2,1,i)*ones(1,numParticle)] + sig_init*randn(2,numParticle);
    weight = ones(1,numParticle)/numParticle;
    for k = 2:1:length
        
        %%% time update
        for n = 1:1:numParticle
            particle(:,n) = particle(:,n) + v_meas(:,k-1,i)*dt;
            particle(:,n) = particle(:,n) + sig_proc*randn(2,1)*dt;
        end
        
        cov_mi = cov(particle');
        dcov_mi = det(cov_mi);
  
        %%% measurement update
        for n = 1:1:numParticle      
            z_est = DEM_height(particle(:,n),DEM);
            weight(n) = weight(n) * likelihood(z_est, z_meas(k,i), sig_meas);
        end
        weight = weight/sum(weight);
        
        %%% resampling (SIR)
        [particle, weight] = Resample(particle, weight);
        
        x_est(:,k,i) = mean(particle'); % MLSE estimate
        x_err(:,k,i) = x_est(:,k,i) - x_true(:,k);
      
        %%% cov test
        dcov_pl = det(cov(particle'));

        dcov_res(:,k,i) = covAnal(cov(particle'));
        
%         [numMode(k), significance, dist(k)] = NumMode(particle);
%         signi(k) = significance(1);
        
        if (dcov_pl > dcov_mi) % if posterior cov > prior cov
            covIncrease(k) = covIncrease(k) + 1;
        end
    end
%     if (mod(i,10) == 0)
        disp([' PF Monte-Carlo step ' int2str(i) '/' int2str(numMonte) ' completed']);
%     end
end
disp(['-- PF completed --']);
x_est_PF = x_est;
x_err_PF = x_err;
particle_PF = particle;
covIncrease_PF = covIncrease;
dcov_res_PF = dcov_res;

[X_err_RMS_PF, Y_err_RMS_PF, d_err_RMS_PF] = RMSE(x_err_PF, length, numMonte);

end

%% MPF

if (RUN_MPF == 1)

x_est = zeros(2,length,numMonte);
x_err = zeros(2,length,numMonte);
particle = zeros(2,numParticle);
dcov_res = zeros(2,length,numMonte);
dcov_mi = zeros(length,numMonte);
dcov_pl = zeros(length,numMonte);
numMode = zeros(1,length);
signi = zeros(1,length);
covIncrease = zeros(1,length);
gradx = zeros(1,length);
grady = zeros(1,length);
disp(['-- MPF started --']);
for i = 1:1:numMonte
    rng('shuffle');
    
    x_est(:,1,i) = x_init + x_init_err;
    x_err(:,1,i) = x_init_err;
    particle = [x_est(1,1,i)*ones(1,numParticle); x_est(2,1,i)*ones(1,numParticle)] + sig_init*randn(2,numParticle);
    weight = ones(1,numParticle)/numParticle;
    for k = 2:1:length
        
        %%% time update
        for n = 1:1:numParticle
            particle(:,n) = particle(:,n) + v_meas(:,k-1,i)*dt;
            particle(:,n) = particle(:,n) + sig_proc*randn(2,1)*dt;
        end
        
        cov_mi = cov(particle');
        dcov_mi = det(cov_mi);
 
%         %%% gradient
%         s = 1;
%         stat = mean(particle')';
%         z_x1 = DEM_height(stat+[-s;0], DEM);
%         z_x2 = DEM_height(stat+[s;0], DEM);
%         z_y1 = DEM_height(stat+[0;-s], DEM);
%         z_y2 = DEM_height(stat+[0;s], DEM);
% 
%         H = [(z_x2 - z_x1)/(2*s), (z_y2 - z_y1)/(2*s)];
%         
%         gradx(k) = abs(H(1));
%         grady(k) = abs(H(2));
    
        %%% clustering
        [numMode(k), significance, center, BW_cr] = NumMode(particle);
        cluster = Cluster(particle, center{numMode(k)});
        signi(k) = significance(1);

        %%% measurement update and resampling
        numCluster = numMode(k);
        weight_c = cell(numCluster, 1);
        numParticle_c = zeros(1,numCluster);
        clusterWeight = zeros(1,numCluster);
        x_est_c = zeros(2,numCluster);
        for c = 1:1:numCluster
            numParticle_c(c) = size(cluster{c},2);
            weight_c{c} = ones(1,numParticle_c(c))/numParticle;
            
            for n = 1:1:numParticle_c(c)
                z_est = DEM_height(cluster{c}(:,n),DEM);
                weight_c{c}(n) = weight_c{c}(n) * likelihood(z_est, z_meas(k,i), sig_meas);
            end
            
            clusterWeight(c) = sum(weight_c{c});
            weight_c{c} = weight_c{c}/clusterWeight(c);
            
            [cluster{c}, tmp] = Resample(cluster{c}, weight_c{c});
            x_est_c(:,c) = mean(cluster{c}');   
        end
        
        %%% merge estimates
        x_est(:,k,i) = zeros(2,1);
        clusterWeight = clusterWeight/sum(clusterWeight);
        for c = 1:1:numCluster
            x_est(:,k,i) = x_est(:,k,i) + x_est_c(:,c)*clusterWeight(c);
        end
        x_err(:,k,i) = x_est(:,k,i) - x_true(:,k);
        
        %%% merge particles
        particle = [];
        for c = 1:1:numCluster
            particle = [particle, cluster{c}];
        end
        size(particle,2);
        weight = ones(1,numParticle)/numParticle;
        
        %%% cov test
        dcov_pl = det(cov(particle'));
        dcov_res(:,k,i) = covAnal(cov(particle'));

        if (dcov_pl > dcov_mi) % if posterior cov > prior cov
            covIncrease(k) = covIncrease(k) + 1;
        end
    end
%     if (mod(i,10) == 0)
        disp([' MPF Monte-Carlo step ' int2str(i) '/' int2str(numMonte) ' completed']);
%     end
end
disp(['-- MPF completed --']);
x_est_MPF = x_est;
x_err_MPF = x_err;
particle_MPF = particle;
covIncrease_MPF = covIncrease;
dcov_res_MPF = dcov_res;

[X_err_RMS_MPF, Y_err_RMS_MPF, d_err_RMS_MPF] = RMSE(x_err_MPF, length, numMonte);

end

%% APF

if (RUN_APF == 1)

x_est = zeros(2,length,numMonte);
x_err = zeros(2,length,numMonte);
particle = zeros(2,numParticle);
dcov_res = zeros(2,length,numMonte);
covIncrease = zeros(1,length);
disp(['-- APF started --']);
for i = 1:1:numMonte
    rng('shuffle');
    
    x_est(:,1,i) = x_init + x_init_err;
    x_err(:,1,i) = x_init_err;
    particle = [x_est(1,1,i)*ones(1,numParticle); x_est(2,1,i)*ones(1,numParticle)] + sig_init*randn(2,numParticle);
    weight = ones(1,numParticle)/numParticle;
    for k = 2:1:length
        
        %%% auxiliary
        for n = 1:1:numParticle
            particle(:,n) = particle(:,n) + v_meas(:,k-1,i)*dt;
            z_est = DEM_height(particle(:,n),DEM);
            weight(n) = weight(n) * likelihood(z_est, z_meas(k,i), sig_meas);
        end
        weight = weight/sum(weight);      
        
        %%% resampling (SIR)
        [particle, weight] = Resample(particle, weight);
       
        %%% time update
        for n = 1:1:numParticle
            particle(:,n) = particle(:,n) + sig_proc*randn(2,1)*dt;
        end
        
        dcov_mi = det(cov(particle'));

        %%% measurement update
        for n = 1:1:numParticle      
            z_est = DEM_height(particle(:,n),DEM);
            weight(n) = weight(n) * likelihood(z_est, z_meas(k,i), sig_meas);
        end   
        weight = weight/sum(weight);
        x_est(:,k,i) = [dot(weight,particle(1,:)); dot(weight,particle(2,:))]; % MLSE estimate
        x_err(:,k,i) = x_est(:,k,i) - x_true(:,k);  
        
        %% cov test
        [particle_t, weight_t] = Resample(particle, weight);
        dcov_pl = det(cov(particle_t'));
        
        dcov_res(:,k,i) = covAnal(cov(particle_t'));
        
        if (dcov_pl > dcov_mi) % if posterior cov > prior cov
            covIncrease(k) = covIncrease(k) + 1;
        end
        
    end
%     if (mod(i,10) == 0)
        disp([' APF Monte-Carlo step ' int2str(i) '/' int2str(numMonte) ' completed']);
%     end
end
disp(['-- APF completed --']);

x_est_APF = x_est;
x_err_APF = x_err;
particle_APF = particle;
covIncrease_APF = covIncrease;
dcov_res_APF = dcov_res;

[X_err_RMS_APF, Y_err_RMS_APF, d_err_RMS_APF] = RMSE(x_err_APF, length, numMonte);

end


%% save
save('result.mat');


%%
figure;
plot(signi, 'k:'); hold on;
axis([10 70 0 1]);
for k = 1:1:70
    if (covIncrease_PF(k) > 0)
        plot(k,0, 'ko', 'linewidth', 1.5);
    end
end
ylabel('Probability');
xlabel('Time (s)');
legend('Prob', 'Ambi');

save('result_mode.mat', 'signi', 'covIncrease_PF');

%% plot
figure;
plot(time, d_err_RMS_PF, 'b'); hold on;
plot(time, d_err_RMS_APF, 'r'); hold on;
plot(time, d_err_RMS_OOSM, 'k'); hold on;
plot(time, d_err_RMS_MPF, 'c'); hold on;
xlabel('Time (s)');
ylabel('RMSE (m)');
legend('PF', 'APF', 'OOSM', 'MPF');
%%
ePF = sum(d_err_RMS_PF(51:151))/100
eAPF = sum(d_err_RMS_APF(51:151))/100
eOOSM = sum(d_err_RMS_OOSM(51:151))/100
eMPF = sum(d_err_RMS_MPF(51:151))/100

figure;
plot(time(2:end), mean(dcov_res_PF(1,2:end,:),3), 'b'); hold on;
plot(time(2:end), mean(dcov_res_APF(1,2:end,:),3), 'r');
plot(time(2:end), mean(dcov_res_OOSM(1,2:end,:),3), 'k'); hold on;
plot(time(2:end), mean(dcov_res_MPF(1,2:end,:),3), 'c'); hold on;
xlabel('Time (s)');
ylabel('COV (m)');
legend('PF', 'APF', 'OOSM', 'MPF');



