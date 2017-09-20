% Select only cortical samples
cd ('data/genes/processedData'); 
load('MicroarrayDataWITHcustVariance360DistThresh2_CoordsAssigned.mat')
doNormalise = true;
doNormalScale = false; 
Lcortex = 1:180;
D = cell(6,1);
subjNr = cell(6,1);

for s=1:6
    data = DataExpression{s};
    subjNr{s} = data(:,2);
    if doNormalise
        expData = BF_NormalizeMatrix(data(:,3:end),'scaledRobustSigmoid');
    else
        expData = data(:,3:end);
    end
    select = ismember(data(:,2), Lcortex);
    D{s} = expData(select==1,:);
    subjNr{s} = data(select==1,1);
 
end

expression = vertcat(D{1}, D{2}, D{3}, D{4}, D{5}, D{6}); 
subjects = vertcat(subjNr{1}, subjNr{2}, subjNr{3}, subjNr{4}, subjNr{5}, subjNr{6}); 


if doNormalScale
expression = 2.^(expression); 
end

[W,score,~,~,explained] = pca(expression); 
x = score(:,1); y = score(:,2); z = score(:,3); 
C = subjects; 
S = ones(length(subjects),1)+30; 

figure; h = scatter(x,y,S,C); %'filled','MarkerEdgeColor',[0 .7 .7]);
set(gcf,'color','w');
xlabel(sprintf('PC1, explains %d%% variance', round(explained(1))));
ylabel(sprintf('PC2, explains %d%% variance', round(explained(2))));

figure; h = scatter3(x,y,z,S,C);
xlabel(sprintf('PC1, explains %d%% variance', round(explained(1))));
ylabel(sprintf('PC2, explains %d%% variance', round(explained(2))));
zlabel(sprintf('PC3, explains %d%% variance', round(explained(3))));
% on normal scaled data
% [W,score,~,~,explained] = pca(expression2); %score=score'; W=W';
% x = score(:,1); y = score(:,2); z = score(:,3); 
% C = Subj; 
% S = ones(length(Subj),1)+20; 
% 
% figure; h = scatter(x,y,S,C);
% xlabel(sprintf('PC1, explains %d%% variance', round(explained(1))));
% ylabel(sprintf('PC2, explains %d%% variance', round(explained(2))));
% 
% figure; h = scatter3(x,y,z,S,C);
% xlabel(sprintf('PC1, explains %d%% variance', round(explained(1))));
% ylabel(sprintf('PC2, explains %d%% variance', round(explained(2))));
% zlabel(sprintf('PC3, explains %d%% variance', round(explained(3))));
% % on per-subject normalised data
% 
% [W,score,~,~,explained] = pca(expressionNorm); %score=score'; W=W';
% x = score(:,1); y = score(:,2); z = score(:,3); 
% C = Subj; 
% S = ones(length(Subj),1)+20; 
% 
% figure; h = scatter(x,y,S,C);
% xlabel(sprintf('PC1, explains %d%% variance', round(explained(1))));
% ylabel(sprintf('PC2, explains %d%% variance', round(explained(2))));
% 
% figure; h = scatter3(x,y,z,S,C);
% xlabel(sprintf('PC1, explains %d%% variance', round(explained(1))));
% ylabel(sprintf('PC2, explains %d%% variance', round(explained(2))));
% ylabel(sprintf('PC2, explains %d%% variance', round(explained(3))));
% 
% 
% % 
% % figure; 
% % for i=1:length(Subj)
% %     if Subj(i) == 1
% %         scatter(score(i,1),score(i,2),'.', 'k');hold on; 
% %     elseif Subj(i) == 2
% %         scatter(score(i,1),score(i,2),'.', 'r');hold on; 
% %     elseif Subj(i) == 3
% %         scatter(score(i,1),score(i,2),'.', 'g');hold on; 
% %     elseif Subj(i) == 4
% %         scatter(score(i,1),score(i,2),'.', 'm');hold on; 
% %     elseif Subj(i) == 5
% %         scatter(score(i,1),score(i,2),'.', 'b'); hold on; 
% %     elseif Subj(i) == 6
% %         scatter(score(i,1),score(i,2),'.', 'c'); hold on; 
% %     end
% % end
% = [W, pc] 