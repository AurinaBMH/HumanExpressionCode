% compare correlation between batch-normalised and non-normallised data 

% 1. correlate separately normalised data using sigmoid and limma
% normalised data + sigmoid
% import limma normalised data
normalisedExpression = importNormalisedExpression('normalisedExpression.txt'); 
normExpression = normalisedExpression'; 
% sigmoid normalise it
expNorm1 = BF_NormalizeMatrix(normExpression,'scaledRobustSigmoid');

% separately normalise expression
doNormalise = true; 
Lcortex = 1:34;
D = cell(6,1);
subjNr = cell(6,1);

for s=1:6
    data = DataExpression{s}; 
    select = ismember(data(:,2), Lcortex);
    cortexData = data(select==1,3:end);
    if doNormalise
        expData = BF_NormalizeMatrix(cortexData,'scaledRobustSigmoid');
    else
        expData = cortexData;
    end

    D{s} = expData; 
    subjNr{s} = data(select==1,1);
    R{s} = data(select==1,2);
end

expNorm2 = vertcat(D{1}, D{2}, D{3}, D{4}, D{5}, D{6}); 

% if ti was nor normalised separately, normalise now all together
if ~doNormalise
    expNorm2 = BF_NormalizeMatrix(expNorm2,'scaledRobustSigmoid');
end
% calculate correlation for each gene and average
for g=1:size(expNorm2,2)
    r(g) = corr(expNorm1(:,g), expNorm2(:,g), 'type', 'Spearman'); 
end
avR = mean(r)
figure; hist(r, 100); 









subjects = vertcat(subjNr{1}, subjNr{2}, subjNr{3}, subjNr{4}, subjNr{5}, subjNr{6}); 

nice_cmap = [make_cmap('steelblue',50,30,0);flipud(make_cmap('orangered',50,30,0))];
p = randperm(length(subjects)); 
[W,score,~,~,explained] = pca(expNorm2, 'NumComponents',4);
x = score(p,1); y = score(p,2); z = score(p,3); 
C = subjects(p); 
S = ones(length(subjects),1)+60;
A = zeros(size(C,1),3); 

for s=1:length(C)
    if C(s)==1
        A(s,:) = [0 0 .5]; % dark blue
    elseif C(s)==2
        A(s,:) = [.53 .83 .97]; % light blue
    elseif C(s)==3
        A(s,:) = [.7 .7 .7]; % grey
    elseif C(s)==4
        A(s,:) = [1 .60 .40]; % orange
    elseif C(s)==5
        A(s,:) = [.81 .07 .15]; % red
    elseif C(s)==6
        A(s,:) = [.99 .87 .09]; % yellow
    end
end

figure; h = scatter(x,y,S,A,'filled','MarkerEdgeColor',[.55 .55 .55],'LineWidth',1.5); 
set(gcf,'color','w');
xlabel(sprintf('PC1, explains %d%% variance', round(explained(1))));
ylabel(sprintf('PC2, explains %d%% variance', round(explained(2))));
set(gca,'fontsize',15)




expData = expression; 
expDatadeBatched = BF_NormalizeMatrix(normExpression,'scaledRobustSigmoid');

figure; scatter(expDatadeBatched(:,1), expData(:,1))
