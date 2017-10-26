load('DSvariance.mat') % - generated using S5 script (probes chosen based on variance)
DSscores{1} = DS; probeSelection{1} = 'variance'; 

load('DSpc.mat') % - generated using S5 script (probes chosen based on PC)
DSscores{2} = DS; probeSelection{2} = 'pc'; 

load('DSnoise.mat')
DSscores{3} = DS; probeSelection{3} = 'noise'; 

load('DSmean.mat');
DSscores{4} = DS; probeSelection{4} = 'mean'; 


sz=10; 
figure; title ('Correlation between DS scores'); 
r = zeros(4,4); 
p = zeros(4,4); 
f=1; 
for i=1:4
    for j=i+1:4

        subplot(2,3,f); scatter(DSscores{i}, DSscores{j}, sz, 'filled');
        [r(i,j),p(i,j)] = corr(DSscores{i}', DSscores{j}', 'type', 'Spearman'); 
        xlabel(sprintf('DS for probes selected based on %s', probeSelection{i})); 
        ylabel(sprintf('DS for probes selected based on %s', probeSelection{j}));
        f=f+1; 
    end
end

