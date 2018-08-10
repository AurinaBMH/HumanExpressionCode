% correlate variance and intensity for each probe - we expect negative
% correlation
close all; 
clear all; 
cd ('data/genes/processedData')
load('MicroarrayDataWITHcustProbesUpdatedXXX.mat')

% calculate
doCV = false; 

if doCV
    filterINT = 2; 
    nameLabel = 'coeficient of variation'; 
else
    filterINT = 3; 
    nameLabel = 'variance'; 
end

for i=1:size(Expressionall,1)
    if doCV
        V(i) = (std(Expressionall(i,:)))/mean(Expressionall(i,:));
    else
        V(i) = var(Expressionall(i,:));
    end
    
    I(i) = mean(Expressionall(i,:));
    
end

figure; scatter(I,V); xlabel('Intensity');ylabel(nameLabel); 

[rho,pval] = corr(I',V', 'type', 'Spearman');
[xThresholds,yMeans] = BF_PlotQuantiles(I,V,250,0,1); 
xlabel('Mean probe intensity');ylabel(sprintf('Probe %s', nameLabel));
set(gca,'FontSize', 18)


% filetr out probes with average intensity lower than 3
lowIND = find(I<filterINT);
lowInt = I(lowIND);

for i=1:length(lowIND)
    if doCV
        Vlow(i) = (std(Expressionall(lowIND(i),:)))/ mean(Expressionall(lowIND(i),:));
    else
        Vlow(i) = var(Expressionall(lowIND(i),:));
    end
    Ilow(i) = mean(Expressionall(lowIND(i),:));
 
end

figure; scatter(Ilow,Vlow); xlabel('Intensity');ylabel(nameLabel); 


[rhoLOW,pvalLOW] = corr(Ilow',Vlow', 'type', 'Spearman');
[xThresholds,yMeans] = BF_PlotQuantiles(Ilow,Vlow,200,0,1); xlabel('Intensity');ylabel(nameLabel); 


highIND = find(I>=filterINT);
highInt = I(highIND);

for i=1:length(highIND)
    
    if doCV
        Vhigh(i) = (std(Expressionall(highIND(i),:)))/ mean(Expressionall(highIND(i),:));
    else
        Vhigh(i) = var(Expressionall(highIND(i),:));
    end

    Ihigh(i) = mean(Expressionall(highIND(i),:));
    
    %     V(i) = var(expressionAll{3}(:,i));
    %     I(i) = mean(expressionAll{3}(:,i));
end
figure; scatter(Ihigh,Vhigh); xlabel('Intensity');ylabel(nameLabel); 

[rhoHIGH,pvalHIGH] = corr(Ihigh',Vhigh', 'type', 'Spearman');
[xThresholds,yMeans] = BF_PlotQuantiles(Ihigh,Vhigh,200,0,1); xlabel('Intensity');ylabel(nameLabel); 
