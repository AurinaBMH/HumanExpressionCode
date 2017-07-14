clear all; close all;

Parcellation = {'cust250'}; 
Threshold = 2; 
NormMethod = {'zscore'};
LEFTcortex = 1; 
% choose 1 if want to normalise samples assigned to left cortex separately; 
% choose 2 if want to normalise LEFT cortex + left subcortex together
% choose 3 if you want to normalise the whole brain. 


if strcmp(Parcellation, 'aparcaseg')
     
                NumNodes = 82;
                LeftCortex = 34;
                LeftSubcortex = 41; 
                RightCortex = 75;
                RightSubcortex = NumNodes;
            elseif strcmp(Parcellation, 'cust100')
                NumNodes = 220;
                LeftCortex = 100;
                LeftSubcortex = 110; 
                RightCortex = 210;
                RightSubcortex = NumNodes;
            elseif strcmp(Parcellation, 'cust250')
                NumNodes = 530;
                LeftCortex = 250;
                LeftSubcortex = 265; 
                RightCortex = 515;
                RightSubcortex = NumNodes;

 end

load(sprintf('%d_DistThresh%d_%s_combined_ExpressionProbePCA_GeneThr0.mat', NumNodes, Threshold, Parcellation{1}));

 ExpressionSubjROI = cell(6,1);
 CoordinatesSubjROI = cell(6,1);
 
for i=1:6
  % normalise data for each subject separately using samples  
    ExpSubj1 = DataExpression{i,1};
    Coord1 = DataCoordinates{i,1};
    
    if LEFTcortex == 1
        ExpSubj = ExpSubj1((ExpSubj1(:,2)<=LeftCortex),:);
        Coord = Coord1((Coord1(:,5)<=LeftCortex),2:4);
    elseif LEFTcortex == 2
        ExpSubj = ExpSubj1((ExpSubj1(:,2)<=LeftSubcortex),:);
        Coord = Coord1((Coord1(:,5)<=LeftSubcortex),2:4);
    elseif LEFTcortex == 3
        ExpSubj = ExpSubj1;
        Coord = Coord1(:,2:4);
   elseif LEFTcortex == 4
        ExpSubj3 = ExpSubj1((ExpSubj1(:,2)>=LeftSubcortex & ExpSubj1(:,2)<=RightCortex),:);
        Coord3 = Coord1((Coord1(:,5)>=LeftSubcortex & Coord1(:,5)<=RightCortex),2:4);
        
        ExpSubj2 = ExpSubj1((ExpSubj1(:,2)<=LeftCortex),:);
        Coord2 = Coord1((Coord1(:,5)<=LeftCortex),2:4);
        
        ExpSubj = cat(1, ExpSubj2, ExpSubj3);
        Coord = cat(1, Coord2, Coord3); 
    end
    data = ExpSubj(:,3:size(ExpSubj,2));
    %data = ExpSubj(:,3:size(DataExpressionCombined,2));
    DataNORM = BF_NormalizeMatrix(data, NormMethod{1}); 
  % average normalised samples within a ROI for each subject separately
    Inf = ExpSubj(:,1:2); 
    DataExpNORM = cat(2, Inf, DataNORM);

    ROIs = linspace(1,length(DataExpNORM(:,2)), length(DataExpNORM(:,2))); 
    ROIs = ROIs'; 
    ExpressionROI = zeros(length(ROIs),size(data,2));
    CoordinatesROI = zeros(length(ROIs),3);
    
%           for j=1:length(ROIs)
%               IndForROI = find(DataExpNORM(:,2)==(ROIs(j)));
%               NumberOfProbes = length(IndForROI);
%                    fprintf(1,'%u samples for %u ROI found \n', NumberOfProbes, ROIs(j))      
%              % take expression values for a selected entrezID
%              ExpressionRepIntensity = DataExpNORM(IndForROI,3:size(ExpSubj,2));
%              CoordinatesRepIntensity = Coord(IndForROI,:);
% 
%              % calculate the mean for expression data for a selected entrezID
%              ExpressionROI(j,:)= ExpressionRepIntensity;
%              CoordinatesROI(j,:)= CoordinatesRepIntensity;
%           end
      S = zeros(length(ROIs),1);
      S(:) = i; 
      ExpressionSubjROI{i} = [ROIs, S,DataNORM];
      CoordinatesSubjROI{i} = [ROIs, S, Coord];
      
end

CombinedExp = cat(1,ExpressionSubjROI{1}, ExpressionSubjROI{2}, ExpressionSubjROI{3}, ExpressionSubjROI{4}, ExpressionSubjROI{5}, ExpressionSubjROI{6});
CombinedCoord = cat(1,CoordinatesSubjROI{1}, CoordinatesSubjROI{2}, CoordinatesSubjROI{3}, CoordinatesSubjROI{4}, CoordinatesSubjROI{5}, CoordinatesSubjROI{6});
%% plot normalised combined data clustered according to provile similarity. 

DATA = CombinedExp(:,3:size(ExpressionROI,2)); 
Dist = pdist2(DATA, DATA);
d_row = pdist(Dist,'corr');
links_row = linkage(d_row);
[~,~,ord_row] = dendrogram(links_row,0); 

ROI = CombinedExp(:,1); 
Subj = CombinedExp(:,2); 

figure; subplot(1,25,1); imagesc(Subj(ord_row)); 
        subplot(1,25,2); imagesc(ROI(ord_row)); 
        subplot(1,25, [3 25]);  imagesc(DATA(ord_row,:)); caxis([0 1]); title(sprintf('LeftCortex %s', NormMethod{1}));

%% average ROIs between subjects (2nd level average)
% ROIs2 = unique(ROI);        
% ExpressionROIall = zeros(length(ROIs2), size(CombinedExp,2)-2);
% CoordinatesROIall = zeros(length(ROIs2), 3);
% 
%         for k=1:length(ROIs2)
%                                 IndForROI2 = find(CombinedExp(:,1)==(ROIs2(k)));
%                                 NumberOfProbes2 = length(IndForROI2);
%                                    fprintf(1,'%u samples for %u ROI \n', NumberOfProbes2, ROIs2(k))      
%                              % take expression values for a selected entrezID
% 
%                              ExpressionRepIntensity2 = CombinedExp(IndForROI2,3:size(CombinedExp,2));
%                              CoordinateRepIntensity2 = CombinedCoord(IndForROI2,3:5);
%                              
%                             ExpressionROIall(k,:)= mean(ExpressionRepIntensity2,1);
%                             CoordinatesROIall(k,:) = mean(CoordinateRepIntensity2,1);
%                             
%         end
%         
%% calculate and plot coexpression VS distance relationships. 

% ExpressionNorm = corrcoef(ExpressionROIall');
% 
% %make a vector for coexpression and distances
%  Rvect = triu(ExpressionNorm,1);
%  Rvect = (Rvect(:));
% 
%  MRIvoxCoordinates = pdist2(CoordinatesROIall, CoordinatesROIall);
%  Dvect = triu(MRIvoxCoordinates,1);
%  Dvect = (Dvect(:));
% 
% %Plot quantiles
%  BF_PlotQuantiles(Dvect,Rvect,100,0,1)