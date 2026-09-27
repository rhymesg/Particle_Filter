function verify_auxiliary
% Check APF importance ratios and selected-ancestor correspondence.
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
old_path = path;
old_rng = rng;
cleanup = onCleanup(@() restore(old_path, old_rng)); %#ok<NASGU>
addpath(root);
rng(7, 'twister');
particles = [ones(1,1000), 2*ones(1,1000)];
lookahead = [0.9*ones(1,1000), 0.1*ones(1,1000)];
[selected, ~, ancestors] = Resample(particles, lookahead);
assert(isequal(selected, particles(ancestors)));
% With no process noise, the second stage must not count data twice.
w = auxiliary_weights(lookahead(ancestors), lookahead, ancestors);
assert(max(abs(w-1/numel(w))) < 1e-14);
assert(abs(sum(w(selected==1))-0.9) < 0.03);
w = auxiliary_weights([0.4, 0.6], [0.2, 0.3], [2,1]);
assert(norm(w-[4/13,9/13], inf) < 1e-12);
[selected, w] = Resample([1,2,3;4,5,6], [0,1,0]);
assert(isequal(selected,repmat([2;5],1,3)));
assert(norm(w-ones(1,3)/3,inf) < 1e-12);
fprintf('Auxiliary particle checks passed.\n');
end
function restore(old_path, old_rng)
path(old_path);
rng(old_rng);
end
