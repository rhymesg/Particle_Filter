function verify_bandwidth
% Check fallback centers and a mode transition on the final search step.
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
original_path = path;
original_folder = pwd;
fixture = tempname;
mkdir(fixture);
cleanup = onCleanup(@() restore(original_path, original_folder, fixture)); %#ok<NASGU>
write_fixture(fullfile(fixture, 'dskensity2d.m'), sprintf([ ...
    'function [pdf, xi, yi] = dskensity2d(varargin)\n' ...
    'global bw_fixture_count bw_fixture_transition\n' ...
    'bw_fixture_count = bw_fixture_count+1;\n' ...
    'xi = [10,20,30]; yi = [30,40,50]; pdf = zeros(3); pdf(1,1)=1;\n' ...
    'if bw_fixture_transition && bw_fixture_count == 101, pdf=eye(3); end\n' ...
    'end\n']));
write_fixture(fullfile(fixture, 'imregionalmax.m'), sprintf( ...
    'function maxima = imregionalmax(pdf)\nmaxima = logical(pdf);\nend\n'));
addpath(root);
cd(fixture);
clear dskensity2d imregionalmax
global bw_fixture_count bw_fixture_transition
particles = [-1,1,-1,1; -1,-1,1,1];
bw_fixture_count = 0;
bw_fixture_transition = false;
[bandwidth, centers] = FindCriticalBW(particles);
assert(norm(bandwidth-[0.1;0.1], inf) < 1e-12);
assert(isequal(centers{1}, [10;30]));
bw_fixture_count = 0;
bw_fixture_transition = true;
[bandwidth, centers] = FindCriticalBW(particles);
expected = 0.1 + (sqrt(4/3)-0.1)/100;
assert(norm(bandwidth-[expected;expected], inf) < 1e-12);
assert(isequal(centers{1}, [10;30]));
fprintf('Critical-bandwidth checks passed.\n');
end

function write_fixture(filename, content)
file = fopen(filename, 'w');
fprintf(file, '%s', content);
fclose(file);
end

function restore(original_path, original_folder, fixture)
cd(original_folder);
path(original_path);
clear dskensity2d imregionalmax
clear global bw_fixture_count bw_fixture_transition
rmdir(fixture, 's');
end
