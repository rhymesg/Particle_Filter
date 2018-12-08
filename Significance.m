function [ significance ] = Significance( particle, BW_cr, m )
%SIGNIFICANCE 
%   assess significance of the critical bandwidth 

iter = 100;
numParticle = size(particle,2);
sig = sqrt(cov(particle'));
sigx = sig(1,1);
sigy = sig(2,2);

x = zeros(1,numParticle);
y = zeros(1,numParticle);
nullhypo = 0;

% if (m == 2)
%     [pdfxy, xi, yi] = dskensity2d(particle, BW_cr);
%     [xxi,yyi]     = meshgrid(xi,yi);
%     figure (1);
% mesh(xxi,yyi,pdfxy)
% set(gca,'XLim',[min(xi) max(xi)])
% set(gca,'YLim',[min(yi) max(yi)])
%     
% end

for i = 1:1:iter;
    
    % % bootstrap sample 
    for n = 1:1:numParticle;
        ind = ceil(rand*numParticle);
        xs = particle(1,ind);
        ys = particle(2,ind);
        
        x(n) = (1 + BW_cr(1)^2/sigx^2)^(-1/2) * (xs + BW_cr(1)*randn);
        y(n) = (1 + BW_cr(2)^2/sigy^2)^(-1/2) * (ys + BW_cr(2)*randn);
    end
    
    [pdfxy, xi, yi] = dskensity2d([x;y], BW_cr);

    % find local maxima
    LM = imregionalmax(pdfxy);
    numMode = sum(sum(LM));
    
    if (numMode > m)
        nullhypo = nullhypo + 1;
    end

end

significance = 1- nullhypo/iter;

end

