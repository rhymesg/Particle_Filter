function [ BW_cr, center ] = FindCriticalBW( particle )
%FINDCRITICALBW 
%   find critical bandwidth for each mode
%   could be improved by using binary search

sig = sqrt(cov(particle'));
BW_max = [sig(1,1); sig(2,2)];
BW_min = [0.1; 0.1];
iter = 100;
numCheckMode = 5; % check up to mode 5

BW = BW_max;

for m = 1:1:numCheckMode

    BW_max = BW;
    BW_prev = BW;
    [pdfxy, xi, yi] = dskensity2d(particle, BW_prev);
    LM = imregionalmax(pdfxy);
    for i = 1:1:iter

        BW_prev = BW;
        LM_prev = LM;
        BW = BW - (BW_max - BW_min)/iter;
        if (BW < 0)
            error('BW below 0');
        end
        
        [pdfxy, xi, yi] = dskensity2d(particle, BW);

        % find local maxima
        LM = imregionalmax(pdfxy);
        numMode = sum(sum(LM));
        
        if (numMode > m)
            BW_cr(:,m) = BW_prev;
            [c, r] = find(LM_prev, m); % x,y ¹Ý´ë·Î
        	center{m} = [xi(r); yi(c)];
            
            if (numMode > m + 1) % finish finding BW_cr
                return;
            end
            
            break;
        end
    end
    if(i == iter)
        BW_cr(:,m) = BW_min;
        center{m} = [xi(r); yi(c)];
        return;
    end
    
end

end

