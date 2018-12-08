function [ mode, significance, center, BW_cr ] = NumMode( particle )
%NUMMODE
%   B. W. Silverman (1981), Using Kernel Density Estimates to Investigate
%   Multimodality - 2d version (yjkim)
%   Input: particle (2 by n matrix)

[BW_cr, center] = FindCriticalBW(particle);

numCheckMode = size(BW_cr,2);

significance = zeros(1,numCheckMode);
for m = 1:1:numCheckMode

    significance(m) = Significance( particle, BW_cr(:,m), m );

end

if (numCheckMode == 1)
    mode = 1;
%     dist = 0;
    return;
end

% % mode as one with the maximum significance
[val, ind] = max(significance);
mode = ind;


end

