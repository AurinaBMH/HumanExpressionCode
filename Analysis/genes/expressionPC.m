% Select only cortical samples
cd ('data/genes/processedData'); 
load('MicroarrayDataWITHcustLessNoise82DistThresh2_CoordsAssigned.mat')
doNormalise = false;
doNormalScale = false; 
Lcortex = 1:34;
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

p = randperm(length(subjects)); 
[W,score,~,~,explained] = pca(expression, 'NumComponents',4);
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

figure; h = scatter3(x,y,z,S,A,'filled','MarkerEdgeColor',[.7 .7 .7],'LineWidth',1.5); 
xlabel(sprintf('PC1, explains %d%% variance', round(explained(1))));
ylabel(sprintf('PC2, explains %d%% variance', round(explained(2))));
zlabel(sprintf('PC3, explains %d%% variance', round(explained(3))));
set(gca,'fontsize',15)
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