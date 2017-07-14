%
clear all; close all;

 
Parcellation = {'cust100'};
Threshold = 2; 
NormMethod = {'zscore'};
LEFTcortex = 1; 
% choose 1 if want to normalise samples assigned to left cortex separately; 
% choose 2 if want to normalise LEFT cortex + left subcortex together
% choose 3 if you want to normalise the whole brain. 
% choose 4 if you want to normalise left cortex + right cortex. 

Thr = 0;

%Fit = {'exp'};
%Choose what proportion of top DS genes to keep
percent = 5;


% choose 1 if want to normalise samples assigned to left cortex separately; choose 0 if want to normalise all samples together. 
if LEFTcortex == 1 || LEFTcortex == 2
    NumSubjects = 6;
elseif LEFTcortex == 3 || LEFTcortex == 4
    NumSubjects = 2;
end

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
cd ('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Microarray/S01-S06_combined/');
load(sprintf('%d_DistThresh%d_%s_combined_ExpressionProbePCA_GeneThr%d.mat', NumNodes, Threshold, Parcellation{1}, round(Thr)));


 ExpressionSubjROI = cell(6,1);
 CoordinatesSubjROI = cell(6,1);
 CoordSample = cell(6,1);
 ExpSampNorm = cell(6,1); 
 
for i=1:NumSubjects
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
    CoordSample{i} = Coord;
    % normalise sample x gene data for each subject separately   
    %% commented - not normalise
      if strcmp(NormMethod, 'hampel')
            ExpressionSampleNorm = Norm_hampel(data);
      else
            ExpressionSampleNorm = BF_NormalizeMatrix(data, NormMethod{1});
      end

    ROI = ExpSubj(:,2); 

    ExpSampNorm{i} = [ROI, ExpressionSampleNorm]; 
    ROIs = unique(ExpSubj(:,2));
  
    
%% do differential stability calculation:
        %1. average samples to ROIs for each subject
        %2. get DS values for each gene
        %3. take top 5% of DS genes
        
%1. average samples to ROIs for each subject

    ExpressionROI = zeros(length(ROIs),size(data,2));
    CoordinatesROI = zeros(length(ROIs),3);
    
          for j=1:length(ROIs)
              IndForROI = find(ExpSubj(:,2)==(ROIs(j)));
              NumberOfProbes = length(IndForROI);
                   fprintf(1,'%u samples for %u ROI found \n', NumberOfProbes, ROIs(j))      
             % take expression values for a selected entrezID
             %% changed into averaging non normalised data ("ExpressionSampleNorm" changed to "data")
             ExpressionRepIntensity = data(IndForROI,:);
             CoordinatesRepIntensity = Coord(IndForROI,:); 

             % calculate the mean for expression data for a selected entrezID
             ExpressionROI(j,:)= mean(ExpressionRepIntensity,1);
             CoordinatesROI(j,:)= mean(CoordinatesRepIntensity,1);
             
          end

       S = zeros(length(ROIs),1);
       S(:) = i; 
       ExpressionSubjROI{i} = [ROIs, S, ExpressionROI];
       CoordinatesSubjROI{i} = [ROIs, S, CoordinatesROI];
end
% combine noramlised data for all subjects. 

ExpSampNormalisedAll = vertcat(ExpSampNorm{1}, ExpSampNorm{2},ExpSampNorm{3},ExpSampNorm{4},ExpSampNorm{5},ExpSampNorm{6});
CombinedCoord = cat(1,CoordSample{1}, CoordSample{2}, CoordSample{3},...
                      CoordSample{4}, CoordSample{5}, CoordSample{6});
%2. get DS values for each gene     
inter = cell(NumSubjects,NumSubjects);
indexj = cell(NumSubjects,NumSubjects);
indexk = cell(NumSubjects,NumSubjects);
indexjp = cell(NumSubjects,NumSubjects);
Corellations = cell(NumSubjects,NumSubjects);
NumGenes = size(ExpressionSubjROI{1,1},2)-2;  
 


%get ROIs that are in all subjects
R = cell(1,NumSubjects);
for o=1:NumSubjects
R{:,o} = ExpressionSubjROI{o}(:,1);
end
Intersect = mintersect(R{1}, R{2}, R{3}, R{4}, R{5}, R{6});

ROIsindex = zeros(length(Intersect),NumSubjects); % use a set list of ROIS that are present in all subjects
for j=1:NumSubjects

    for w=1:length(Intersect)
        ROIsindex(w,j) = find(R{j}==Intersect(w));
    end
end
        

for j=1:NumSubjects
    for k=j+1:NumSubjects
        % use a set list of ROIS that are present in all subjects
%         [inter{j,k}, indexj{j,k}, indexk{j,k}] = intersect(ExpressionSubjROI{j}(:,1), ExpressionSubjROI{k}(:,1));
%         Exp1 = ExpressionSubjROI{j}(indexj{j,k},3:NumGenes+2);
%         Exp2 = ExpressionSubjROI{k}(indexk{j,k},3:NumGenes+2);
        Exp1 = ExpressionSubjROI{j}(ROIsindex(:,j),3:NumGenes+2);
        Exp2 = ExpressionSubjROI{k}(ROIsindex(:,k),3:NumGenes+2);
        Genes = zeros(1,NumGenes);
        for g=1:NumGenes
           
            Genes(g) = corr(Exp1(:,g),Exp2(:,g),'type','Spearman'); 
        end
        Corellations{j,k} = Genes; 
    end
end

if NumSubjects == 2
    C = vertcat(Corellations{1,2});
elseif NumSubjects == 6
C = vertcat(Corellations{1,2}, Corellations{1,3}, Corellations{1,4}, Corellations{1,5}, Corellations{1,6}, ... 
    Corellations{2,3}, Corellations{2,4}, Corellations{2,5}, Corellations{2,6}, Corellations{3,4}, Corellations{3,5}, ... 
    Corellations{3,6}, Corellations{4,5}, Corellations{4,6}, Corellations{5,6});
end

DS = mean(C,1); 

% combined
% CombinedExp = cat(1,ExpressionSubjROI{1},  ExpressionSubjROI{2}, ExpressionSubjROI{3},...
%     ExpressionSubjROI{4}, ExpressionSubjROI{5}, ExpressionSubjROI{6});

%3. take top % of DS genes 

NrGenes = round(length(DS)*percent/100);

    [ b, ix ] = sort( DS(:), 'descend' );

    DSvalues = zeros(NrGenes, 2);
    for ii=1:NrGenes
        DSvalues(ii,2) = b(ii);
        DSvalues(ii,1) = ix(ii);
    end

% get probeIDs for selected DS genes
Probes = ProbeInformation.ProbeName(DSvalues(:,1));
DSProbeTable = table(Probes, DSvalues(:,2));

%% select selected genes and calculate sample - sample coexpression

SelectedGenes = ExpSampNormalisedAll(:,2:end);
% take genes with highest DS values
SelectedGenes = SelectedGenes(:,DSvalues(:,1));
% calculate sample-sample coexpression
SampleCoexpression = corr(SelectedGenes', 'type', 'Spearman'); 
%% check coexpression - distance relationship. 

MRIvoxCoordinates = pdist2(CombinedCoord, CombinedCoord);

%make a vector for coexpression and distances

DistExpVect(:,1) = MRIvoxCoordinates(:);
%b replace diagonal with zeros
%SampleCoexpression(logical(eye(size(SampleCoexpression)))) = 0;
DistExpVect(:,2) = SampleCoexpression(:);
DistExpVect( ~any(DistExpVect,2), : ) = [];  %rows

Dvect = DistExpVect(:,1);
Rvect = DistExpVect(:,2);

figure; imagesc(SampleCoexpression); caxis([-1,1]);title('Sample-sample coexpression');
colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);

% fit distance correction according to a defined rule
%[f_handle,Stats,c] = GiveMeFit(DistExpVect(:,1),DistExpVect(:,2),Fit{1});
% fit sigmoid
[param,stat] = sigm_fit(DistExpVect(:,1),DistExpVect(:,2));

% plot original doexpression-distance .
BF_PlotQuantiles(DistExpVect(:,1),DistExpVect(:,2),50,0,1); title('Coexpresion vs distance'); ylim([-0.8 1]); 
%hold on; plot(c); 
hold on; scatter(DistExpVect(:,1),stat.ypred,1, '.', 'r');


% switch Fit{1}
% 
%         case 'linear'
%             FitCurve = c.p1*Dvect + c.p2;
%         case 'exp'
%             FitCurve = c.A*exp(-c.n*Dvect) + c.B;
%         case 'exp_1_0'
%             FitCurve = exp(-c.n*Dvect);
%         case 'decay'
%             FitCurve = c.A/Dvect + c.B;
%             Residuals = Rvect' - FitCurve;
%         case 'exp0'
%             FitCurve = c.A.*exp(-c.n*Dvect);
%         case 'exp1'
%             FitCurve = exp(-c.n*Dvect) + c.B;
% end
% get residuals
%Residuals = Rvect - FitCurve;
Residuals = Rvect - stat.ypred;
BF_PlotQuantiles(DistExpVect(:,1),nonzeros(Residuals(:)),50,0,1); title('Coexpresion vs distance corrected'); ylim([-0.8 1]); 


%% Plot corrected sample - sample coexpression matrix; 


NumSamples = size(MRIvoxCoordinates,1);
CorrectedCoexpression = reshape(Residuals,[NumSamples, NumSamples]);
figure; imagesc(CorrectedCoexpression); caxis([-1,1]);title('Corrected Sample-sample coexpression');
colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);

%% average coexpression values within a ROI and plot the corrected matirx (ROI-ROI coexpression);
W = unique(ExpSampNormalisedAll(:,1));

ParcelCoexpression = zeros(length(W),length(W));
ROIs = ExpSampNormalisedAll(:,1);

[sROIs, ind] = sort(ROIs);
CorrectedCoexpressionSorted = CorrectedCoexpression(ind, ind);
CoexpressionSorted = SampleCoexpression(ind, ind);


figure; subplot(1,25,[1 2]); imagesc(sROIs);
subplot(1,25,[3 25]); imagesc(CoexpressionSorted); caxis([-1,1]); title('Corrected coexpression sorted samples');
colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);

for i=1:length(W)
    for j=1:length(W)

        A = find(sROIs == W(i));
        B = find(sROIs == W(j));
        %for corrected
        P = CorrectedCoexpressionSorted(A, B);
        %for uncorrected
        %P = CoexpressionSorted(A, B);
        ParcelCoexpression(i,j) = mean(mean(P));

    end
end

figure; imagesc(ParcelCoexpression); caxis([-1,1]); title('Parcellation coexpression ROIs');
colormap([flipud(BF_getcmap('blues',9));[1 1 1];BF_getcmap('reds',9)]);

% add zeros values to missing ROIs; p - missing ROI
p = setdiff(linspace(1,LeftCortex,LeftCortex), unique(ExpSampNormalisedAll(:,1)));  
Exp = ParcelCoexpression;
N1 = nan(LeftCortex,1);
N2 = nan(1, LeftCortex-1);


B = vertcat(Exp(1:p-1,:), N2, Exp(p:end,:));
Exp = horzcat(B(:,1:p-1), N1, B(:,p:end));





