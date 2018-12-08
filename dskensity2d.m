function [ pdfxy, xi, yi ] = dskensity2d( particle, BW )
%DSKENSITY2D 
%   2d kernel density estimate

x = particle(1,:)';
y = particle(2,:)';

% Estimate a continuous pdf from the discrete data
[pdfx, xi]= ksdensity(x, 'bandwidth', BW(1));
[pdfy, yi]= ksdensity(y, 'bandwidth', BW(2));

% Create 2-d grid of coordinates and function values, suitable for 3-d plotting
% [xxi, yyi]     = meshgrid(xi,yi);
[pdfxx, pdfyy] = meshgrid(pdfx,pdfy);
% Calculate combined pdf, under assumption of independence
pdfxy = pdfxx.*pdfyy; 


end

