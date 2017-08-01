%% Author: Aurina

%Last modiffied: 2017-07-31
%Last modiffied: 2017-08-01

%------------------------------------------------------------------------------
% Choose options
%------------------------------------------------------------------------------
useCUSTprobes = false; % choose if you want to use data with CUST probes
probeSelection = 'Variance';% (Variance', LessNoise', 'Mean')
parcellation = 'aparcaseg';%, 'cust100', 'cust250'};
distanceThreshold = 2; % first run 30, then with the final threshold 2
percentDS = 5;
coexpressionFor = 'all'; 
Fit = {'exp'}; 
normMethod = 'scaledRobustSigmoid';
normaliseWhat = 'Lcortex'; %(LcortexSubcortex, wholeBrain, LRcortex)
% choose Lcortex if want to normalise samples assigned to left cortex separately; 
% choose LcortexSubcortex if want to normalise LEFT cortex + left subcortex together
% choose wholeBrain if you want to normalise the whole brain. 
% choose LRcortex if you want to normalise left cortex + right cortex. 

%------------------------------------------------------------------------------
% Define number of subjects and parcellation details based on choises
%------------------------------------------------------------------------------

switch normaliseWhat
    case 'Lcortex'
        subjects = 1:6;
    case 'LcortexSubcortex'
        subjects = 1:6;
    case 'wholeBrain'
        subjects = 1:2;
    case 'LRcortex'
        subjects = 1:2;
end

if strcmp(parcellation, 'aparcaseg')
     
                NumNodes = 82;
                LeftCortex = 34;
                LeftSubcortex = 41; 
                RightCortex = 75;
                RightSubcortex = NumNodes;
                
    elseif strcmp(parcellation, 'cust100')
                NumNodes = 220;
                LeftCortex = 100;
                LeftSubcortex = 110; 
                RightCortex = 210;
                RightSubcortex = NumNodes;
                
    elseif strcmp(parcellation, 'cust250')
                NumNodes = 530;
                LeftCortex = 250;
                LeftSubcortex = 265; 
                RightCortex = 515;
                RightSubcortex = NumNodes;
                
end

cd ('data/genes/processedData');
load(sprintf('MicroarrayDatad%s%dDistThresh%d_CoordsAssigned.mat', probeSelection, NumNodes, distanceThreshold)); 


 expressionSubjROI = cell(6,1);
 coordinatesSubjROI = cell(6,1);
 coordSample = cell(6,1);
 expSampNorm = cell(6,1); 
 
%----------------------------------------------------------------------------------
% Normalise data for each subject separately 
% Do differential stability calculation:
        %1. average expression values for each ROI for differential stability calculation
        %2. get DS values for each gene
        %3. take top 5% of DS genes
%----------------------------------------------------------------------------------

for sub=subjects
  % normalise data for each subject separately using samples  
    expSingleSubj = DataExpression{sub,1};
    coordSingle = DataCoordinates{sub,1};
    
    switch normaliseWhat
        
        case 'Lcortex'
            
        expSubj = expSingleSubj((expSingleSubj(:,2)<=LeftCortex),:);
        coord = coordSingle((coordSingle(:,2)<=LeftCortex),3:5);
        
        case 'LcortexSubcortex'
            
        expSubj = expSingleSubj((expSingleSubj(:,2)<=LeftSubcortex),:);
        coord = coordSingle((coordSingle(:,2)<=LeftSubcortex),3:5);
        
        case 'wholeBrain'
            
        expSubj = expSingleSubj;
        coord = coordSingle(:,3:5);
        
        case 'LRcortex'
            
        expSubjRight = expSingleSubj((expSingleSubj(:,2)>=LeftSubcortex & expSingleSubj(:,2)<=RightCortex),:);
        coordRight = coordSingle((coordSingle(:,2)>=LeftSubcortex & coordSingle(:,2)<=RightCortex),3:5);
        
        expSubjLeft = expSingleSubj((expSingleSubj(:,2)<=LeftCortex),:);
        coordLeft = coordSingle((coordSingle(:,2)<=LeftCortex),3:5);
        
        expSubj = cat(1, expSubjLeft, expSubjRight);
        coord = cat(1, coordLeft, coordRight); 
        
    end
    
    data = expSubj(:,3:size(expSubj,2));
    coordSample{sub} = coord;
    ROI = expSubj(:,2);
    % normalise sample x gene data for each subject separately   
    %% commented - not normalise
    switch normMethod
        case 'hampel'
            dataNorm = Norm_hampel(data);
            fprintf('Normalising gene expression data\n')
        otherwise
            dataNorm = BF_NormalizeMatrix(data, normMethod);
            fprintf('Normalising gene expression data\n')
    end
    
    expSampNorm{sub} = [ROI, dataNorm]; 
    ROIs = unique(expSubj(:,2));

% average expression values for each ROI for differential stability calculation

    expressionROI = zeros(length(ROIs),size(data,2));
    coordinatesROI = zeros(length(ROIs),3);
    
          for j=1:length(ROIs)
              indROI = find(expSubj(:,2)==(ROIs(j)));
              noProbes = length(indROI);
                   fprintf(1,'%u samples for %u ROI found \n', noProbes, ROIs(j))      
             % take expression values for a selected entrezID
             expressionRepInt = data(indROI,:);
             coordinatesRepInt = coord(indROI,:); 

             % calculate the mean for expression data for a selected entrezID
             expressionROI(j,:)= mean(expressionRepInt,1);
             coordinatesROI(j,:)= mean(coordinatesRepInt,1);
             
          end

       S = zeros(length(ROIs),1);
       S(:) = sub; 
       expressionSubjROI{sub} = [S, ROIs, expressionROI];
       coordinatesSubjROI{sub} = [S, ROIs, coordinatesROI];
end
%----------------------------------------------------------------------------------
% Combine noramlised data for all subjects.
%----------------------------------------------------------------------------------

expSampNormalisedAll = vertcat(expSampNorm{1}, expSampNorm{2},expSampNorm{3},expSampNorm{4},expSampNorm{5},expSampNorm{6});
combinedCoord = cat(1,coordSample{1}, coordSample{2}, coordSample{3},...
                      coordSample{4}, coordSample{5}, coordSample{6});

%----------------------------------------------------------------------------------
% Pre-define variables for DS calculation 
%----------------------------------------------------------------------------------

numSubjects = max(subjects); 
inter = cell(numSubjects,numSubjects);
indexj = cell(numSubjects,numSubjects);
indexk = cell(numSubjects,numSubjects);
indexjp = cell(numSubjects,numSubjects);
corellations = cell(numSubjects,numSubjects);
numGenes = size(expressionSubjROI{1,1},2)-2;  

%----------------------------------------------------------------------------------
% Get ROIs that are in all subjects 
%----------------------------------------------------------------------------------
fprintf('Selecting ROIs for DS calculation\n')
R = cell(1,numSubjects);
for o=1:numSubjects
R{:,o} = expressionSubjROI{o}(:,2);
end
intersectROIs = mintersect(R{1}, R{2}, R{3}, R{4}, R{5}, R{6});

% use a set list of ROIS that are present in all subjects
ROIsindex = zeros(length(intersectROIs),numSubjects); 
for j=1:numSubjects

    for w=1:length(intersectROIs)
        ROIsindex(w,j) = find(R{j}==intersectROIs(w));
    end
    
end

%----------------------------------------------------------------------------------
% For each pair of cubjects, calculate correlation (Spearman) of regional expression
% for each gene to select genes that have consistent expression patterns
% through regions between pairs of subjects. 
%----------------------------------------------------------------------------------
fprintf('Calculating differential stability\n')
for j=1:numSubjects
    for k=j+1:numSubjects
        
        expSubone = expressionSubjROI{j}(ROIsindex(:,j),3:numGenes+2);
        expSubtwo = expressionSubjROI{k}(ROIsindex(:,k),3:numGenes+2);
        genes = zeros(1,numGenes);
        
        for g=1:numGenes
            genes(g) = corr(expSubone(:,g),expSubtwo(:,g),'type','Spearman'); 
        end
        
        corellations{j,k} = genes; 
    end
end

%----------------------------------------------------------------------------------
% Combina data for all pairs of subjects
%----------------------------------------------------------------------------------
if numSubjects == 2
    C = vertcat(corellations{1,2});
elseif numSubjects == 6
C = vertcat(corellations{1,2}, corellations{1,3}, corellations{1,4}, corellations{1,5}, corellations{1,6}, ... 
    corellations{2,3}, corellations{2,4}, corellations{2,5}, corellations{2,6}, corellations{3,4}, corellations{3,5}, ... 
    corellations{3,6}, corellations{4,5}, corellations{4,6}, corellations{5,6});
end
%----------------------------------------------------------------------------------
% Take the mean for each gene - this is DS score for a gene
% gene that have most consistent expression pattern through regions will
% get highest scores
%----------------------------------------------------------------------------------
DS = mean(C,1); 

%----------------------------------------------------------------------------------
% Take top % of DS genes 
%----------------------------------------------------------------------------------
fprintf('Selecting genes with highest differential stability \n')
nrGenes = round(length(DS)*percentDS/100);

    [ b, ix ] = sort( DS(:), 'descend' );

    DSvalues = zeros(nrGenes, 2);
    for ii=1:nrGenes
        DSvalues(ii,2) = b(ii);
        DSvalues(ii,1) = ix(ii);
    end
    
%----------------------------------------------------------------------------------
% Get probeIDs for selected DS genes
%----------------------------------------------------------------------------------

probes = probeInformation.ProbeName(DSvalues(:,1));
DSProbeTable = table(probes, DSvalues(:,2));

switch coexpressionFor
    case 'all'
%----------------------------------------------------------------------------------
% Take selected genes and calculate sample - sample coexpression
%----------------------------------------------------------------------------------
selectedGenes = expSampNormalisedAll(:,2:end);
selectedGenes = selectedGenes(:,DSvalues(:,1)); % take genes with highest DS values
SampleCoexpression = corr(selectedGenes', 'type', 'Spearman'); % calculate sample-sample coexpression

%----------------------------------------------------------------------------------
% Check coexpression - distance relationship. 
%----------------------------------------------------------------------------------

MRIvoxCoordinates = pdist2(combinedCoord, combinedCoord);
distExpVect(:,1) = MRIvoxCoordinates(:); % make a vector for distances
SampleCoexpression(logical(eye(size(SampleCoexpression)))) = 0; % replace diagonal with zeros
distExpVect(:,2) = SampleCoexpression(:); % make a vector for coexpression values
distExpVect( ~any(distExpVect,2), : ) = [];  % remove rows for diagonal elememns as thay will have 0 distance and 1 coexpression

Dvect = distExpVect(:,1);
Rvect = distExpVect(:,2);

figure; imagesc(SampleCoexpression); caxis([-1,1]);title('Sample-sample coexpression');
colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);

%----------------------------------------------------------------------------------
% Fit distance correction according to a defined rule
%----------------------------------------------------------------------------------
% 
[f_handle,Stats,c] = GiveMeFit(distExpVect(:,1),distExpVect(:,2),Fit{1});

%[param,stat] = sigm_fit(distExpVect(:,1),distExpVect(:,2));

% plot original coexpression-distance .
BF_PlotQuantiles(distExpVect(:,1),distExpVect(:,2),50,0,1); title('Coexpresion vs distance'); ylim([-0.8 1]); 
%hold on; plot(c); 



switch Fit{1}

        case 'linear'
            FitCurve = c.p1*Dvect + c.p2;
        case 'exp'
            FitCurve = c.A*exp(-c.n*Dvect) + c.B;
        case 'exp_1_0'
            FitCurve = exp(-c.n*Dvect);
        case 'decay'
            FitCurve = c.A/Dvect + c.B;
            Residuals = Rvect' - FitCurve;
        case 'exp0'
            FitCurve = c.A.*exp(-c.n*Dvect);
        case 'exp1'
            FitCurve = exp(-c.n*Dvect) + c.B;
            
end
hold on; scatter(distExpVect(:,1),FitCurve,1, '.', 'r');
% get residuals
Residuals = Rvect - FitCurve;
%Residuals = Rvect - stat.ypred;
BF_PlotQuantiles(distExpVect(:,1),nonzeros(Residuals(:)),50,0,1); title('Coexpresion vs distance corrected'); ylim([-0.8 1]); 



%----------------------------------------------------------------------------------
% Plot corrected sample - sample coexpression matrix
%----------------------------------------------------------------------------------

numSamples = size(MRIvoxCoordinates,1);
CorrectedCoexpression = reshape(Residuals,[numSamples, numSamples]);
figure; imagesc(CorrectedCoexpression); caxis([-1,1]);title('Corrected Sample-sample coexpression');
colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);

%----------------------------------------------------------------------------------
% Plot corrected ROI-ROI coexpression matrix
%----------------------------------------------------------------------------------

W = unique(expSampNormalisedAll(:,1));

ParcelCoexpression = zeros(length(W),length(W));
ROIs = expSampNormalisedAll(:,1);

[sROIs, ind] = sort(ROIs);
CorrectedCoexpressionSorted = CorrectedCoexpression(ind, ind);
CoexpressionSorted = SampleCoexpression(ind, ind);

figure; subplot(1,25,[1 2]); imagesc(sROIs);
subplot(1,25,[3 25]); imagesc(CoexpressionSorted); caxis([-1,1]); title('Corrected coexpression sorted samples');
colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);

for sub=1:length(W)
    for j=1:length(W)

        A = find(sROIs == W(sub));
        B = find(sROIs == W(j));
        %for corrected
        P = CorrectedCoexpressionSorted(A, B);
        %for uncorrected
        %P = CoexpressionSorted(A, B);
        ParcelCoexpression(sub,j) = mean(mean(P));

    end
end

figure; imagesc(ParcelCoexpression); caxis([-1,1]); title('Parcellation coexpression ROIs');
colormap([flipud(BF_getcmap('blues',9));[1 1 1];BF_getcmap('reds',9)]);

% add zeros values to missing ROIs; p - missing ROI
p = setdiff(linspace(1,LeftCortex,LeftCortex), unique(expSampNormalisedAll(:,1)));  
Exp = ParcelCoexpression;
N1 = nan(LeftCortex,1);
N2 = nan(1, LeftCortex-1);


B = vertcat(Exp(1:p-1,:), N2, Exp(p:end,:));
Exp = horzcat(B(:,1:p-1), N1, B(:,p:end));
    case 'separate'
        
        
        
end





