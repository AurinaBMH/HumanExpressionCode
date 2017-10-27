clear all; 
load('DSnewVariance.mat') % - generated using S5 script (probes chosen based on variance)
[a1,order] = sort(probeInformation.EntrezID); 
DSscores{1} = DS(order); probeSelection{1} = 'variance'; 
allData{1} = expSampNormalisedAll(:,order); 

load('DSnewPC.mat') % - generated using S5 script (probes chosen based on PC)
[a2,order] = sort(probeInformation.EntrezID); 
DSscores{2} = DS(order); probeSelection{2} = 'PC'; 
allData{2} = expSampNormalisedAll(:,order); 


load('DSnewLessNoise.mat')
[a3,order] = sort(probeInformation.EntrezID); 
DSscores{3} = DS(order); probeSelection{3} = 'noise'; 
allData{3} = expSampNormalisedAll(:,order); 


load('DSnewMean.mat');
[a4,order] = sort(probeInformation.EntrezID); 
DSscores{4} = DS(order); probeSelection{4} = 'Mean'; 
allData{4} = expSampNormalisedAll(:,order); 

load('DSnewrandom.mat');
[a5,order] = sort(probeInformation.EntrezID); 
DSscores{5} = DS(order); probeSelection{5} = 'random'; 
allData{5} = expSampNormalisedAll(:,order); 


sz=10; 
figure; title ('Correlation between DS scores'); 
r = zeros(5,5); 
p = zeros(5,5); 
f=1; 
for i=1:5
    for j=i+1:5

        subplot(2,5,f); scatter(DSscores{i}, DSscores{j}, sz, 'filled');
        [r(i,j),p(i,j)] = corr(DSscores{i}', DSscores{j}', 'type', 'Spearman'); 
        xlabel(sprintf('DS for probes selected based on %s', probeSelection{i})); 
        ylabel(sprintf('DS for probes selected based on %s', probeSelection{j}));
        f=f+1; 
    end
end

