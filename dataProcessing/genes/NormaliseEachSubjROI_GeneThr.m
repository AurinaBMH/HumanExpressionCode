%clear all; close all;

 
Parcellation = {'aparcaseg'};
Threshold = 2; 
NormMethod = {'scaledRobustSigmoid'};
LEFTcortex = 1; 
Thr = 0;
% choose 1 if want to normalise samples assigned to left cortex separately; 
% choose 2 if want to normalise LEFT cortex + left subcortex together
% choose 3 if you want to normalise the whole brain. 
% choose 4 if you want to normalise left cortex + right cortex. 

% choose 1 if want to normalise samples assigned to left cortex separately; choose 0 if want to normalise all samples together. 

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
cd ('/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S01-S06_combined/');
load(sprintf('%d_DistThresh%d_%s_combined_ExpressionProbePCA_GeneThr%d.mat', NumNodes, Threshold, Parcellation{1}, round(Thr)));


 ExpressionSubjROI = cell(6,1);
 CoordinatesSubjROI = cell(6,1);
for i=1:6
  % normalise data for each subject separately using ROIs
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

     
  % average samples within a ROI for each cubject separately

    ROIs = unique(ExpSubj(:,2)); 
    ExpressionROI = zeros(length(ROIs),size(data,2));
    CoordinatesROI = zeros(length(ROIs),3);
    
          for j=1:length(ROIs)
              IndForROI = find(ExpSubj(:,2)==(ROIs(j)));
              NumberOfProbes = length(IndForROI);
                   fprintf(1,'%u samples for %u ROI found \n', NumberOfProbes, ROIs(j))      
             % take expression values for a selected entrezID
             ExpressionRepIntensity = data(IndForROI,:);
             CoordinatesRepIntensity = Coord(IndForROI,:); 

             % calculate the mean for expression data for a selected entrezID
             ExpressionROI(j,:)= mean(ExpressionRepIntensity,1);
             CoordinatesROI(j,:)= mean(CoordinatesRepIntensity,1);
             
          end
      % normalise ROI x gene data for each subject separately    
      ExpressionROInorm = BF_NormalizeMatrix(ExpressionROI, NormMethod{1});
      S = zeros(length(ROIs),1);
      S(:) = i; 
      ExpressionSubjROI{i} = [ROIs, S, ExpressionROInorm];
      CoordinatesSubjROI{i} = [ROIs, S, CoordinatesROI];
end
% combine averaged, normalised data from all subjects
CombinedExp = cat(1,ExpressionSubjROI{1},  ExpressionSubjROI{2}, ExpressionSubjROI{3}, ExpressionSubjROI{4}, ExpressionSubjROI{5}, ExpressionSubjROI{6});
CombinedCoord = cat(1,CoordinatesSubjROI{1}, CoordinatesSubjROI{2}, CoordinatesSubjROI{3}, CoordinatesSubjROI{4}, CoordinatesSubjROI{5}, CoordinatesSubjROI{6});

%% plot normalised combined data clustered according to provile similarity. 

DATA = CombinedExp(:,3:size(ExpSubj,2)); 
Dist = pdist2(DATA, DATA);
d_row = pdist(Dist,'corr');
links_row = linkage(d_row);
[~,~,ord_row] = dendrogram(links_row,0); 

ROI = CombinedExp(:,1); 
Subj = CombinedExp(:,2); 

figure; subplot(1,25,1); imagesc(Subj(ord_row)); 
        subplot(1,25,2); imagesc(ROI(ord_row)); 
        subplot(1,25, [3 25]);  imagesc(DATA(ord_row,:)); caxis([0 1]); title(sprintf('LeftCortex %s', NormMethod{1}));

        %% average data within ROI between subjects
ROIs2 = unique(ROI);        
ExpressionROIall = zeros(length(ROIs2), size(CombinedExp,2)-2);
CoordinatesROIall = zeros(length(ROIs2), 3);

        for k=1:length(ROIs2)
                                IndForROI2 = find(CombinedExp(:,1)==(ROIs2(k)));
                                NumberOfProbes2 = length(IndForROI2);
                                   fprintf(1,'%u samples for %u ROI \n', NumberOfProbes2, ROIs2(k))      
                             % take expression values for a selected entrezID

                             ExpressionRepIntensity2 = CombinedExp(IndForROI2,3:size(CombinedExp,2));
                             CoordinateRepIntensity2 = CombinedCoord(IndForROI2,3:5);
                             
                            ExpressionROIall(k,:)= mean(ExpressionRepIntensity2,1);
                            CoordinatesROIall(k,:) = mean(CoordinateRepIntensity2,1);
                            
        end
        

ExpressionNorm = corr(ExpressionROIall');
MRIvoxCoordinates = pdist2(CoordinatesROIall, CoordinatesROIall);
 Dvect = triu(MRIvoxCoordinates,1);
 Rvect = triu(ExpressionNorm,1);
%  L = logical(Dvect);
%  Rvect = Rvect.*L;
 Dvect = nonzeros(Dvect(:));
 Rvect = nonzeros(Rvect(:));
 
%make a vector for coexpression and distances


%Plot quantiles
 BF_PlotQuantiles(Dvect,Rvect,100,0,1)
 figure; imagesc(ExpressionNorm); caxis([-1,1])
colormap([flipud(BF_getcmap('blues',9));BF_getcmap('reds',9)]);
        

