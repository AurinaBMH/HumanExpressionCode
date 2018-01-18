% correlate variance and intensity for each probe - we expect negative
% correlation

% calculate 

for i=1:size(Expressionall,1)
    V(i) = var(Expressionall(i,:)); 
    I(i) = mean(Expressionall(i,:)); 
    
%     V(i) = var(expressionAll{3}(:,i)); 
%     I(i) = mean(expressionAll{3}(:,i)); 
end

figure; scatter(I,V); xlabel('Intensity'); ylabel('Variance');
[rho,pval] = corr(I',V', 'type', 'Spearman'); 
[xThresholds,yMeans] = BF_PlotQuantiles(I,V,100,0,1); xlabel('Intensity'); ylabel('Variance');

% filetr out probes with average intensity lower than 3
lowIND = find(I<3); 
lowInt = I(lowIND); 

for i=1:length(lowIND)
    Vlow(i) = var(Expressionall(lowIND(i),:)); 
    Ilow(i) = mean(Expressionall(lowIND(i),:)); 
    
%     V(i) = var(expressionAll{3}(:,i)); 
%     I(i) = mean(expressionAll{3}(:,i)); 
end

figure; scatter(Ilow,Vlow); xlabel('Intensity'); ylabel('Variance');
[rhoLOW,pvalLOW] = corr(Ilow',Vlow', 'type', 'Spearman'); 
[xThresholds,yMeans] = BF_PlotQuantiles(Ilow,Vlow,200,0,1); xlabel('Intensity'); ylabel('Variance');


highIND = find(I>=3); 
highInt = I(highIND); 

for i=1:length(highIND)
    Vhigh(i) = var(Expressionall(highIND(i),:)); 
    Ihigh(i) = mean(Expressionall(highIND(i),:)); 
    
%     V(i) = var(expressionAll{3}(:,i)); 
%     I(i) = mean(expressionAll{3}(:,i)); 
end

figure; scatter(Ihigh,Vhigh); xlabel('Intensity'); ylabel('Variance');
[rhoHIGH,pvalHIGH] = corr(Ihigh',Vhigh', 'type', 'Spearman'); 
[xThresholds,yMeans] = BF_PlotQuantiles(Ihigh,Vhigh,200,0,1); xlabel('Intensity'); ylabel('Variance');