function [expPlot, parcelCoexpression, Residuals, distExpVect, distPlot] = calculateCoexpression(sampleDistances, selectedGenes, DSvalues, W, ROIs,nROIs, Fit, correctDistance, resolution)

%sampleDistances = maskHalf(sampleDistances);

        selectedGenesN = selectedGenes(:,DSvalues(:,1)); % take genes with highest DS values
        sampleCoexpression = corr(selectedGenesN', 'type', 'Spearman'); % calculate sample-sample coexpression
        %sampleCoexpression = maskHalf(sampleCoexpression);
        
        sampleCoexpression(logical(eye(size(sampleCoexpression)))) = NaN; % replace diagonal with NaN

switch resolution
    case 'sample'
        distExpVect(:,1) = sampleDistances(:); % make a vector for distances
        distExpVect(:,2) = sampleCoexpression(:); % make a vector for coexpression values
        distExpVect(any(isnan(distExpVect), 2), :) = [];
        
        %distExpVect(isnan(distExpVect,2), : ) = [];  % remove rows for diagonal elememns as thay will have 0 distance and 1 coexpression
        
        Dvect = distExpVect(:,1);
        Rvect = distExpVect(:,2);
        
        figure; imagesc(sampleCoexpression); caxis([-1,1]);title('Sample-sample coexpression');
        colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);
    case 'ROI'
        [sROIs, ind] = sort(ROIs);
        coexpressionSorted = sampleCoexpression(ind, ind);
        distancesSorted = sampleDistances(ind, ind); 
        
        parcelCoexpression = zeros(length(W),length(W));
        parcelDistances = zeros(length(W),length(W));


        for sub=1:length(W)
            for j=1:length(W)
                
                A = sROIs == W(sub);
                B = sROIs == W(j);

               P = coexpressionSorted(A, B);

                D = distancesSorted(A,B);
                parcelDistances(sub,j) = nanmean(mean(D));
                parcelCoexpression(sub,j) = nanmean(mean(P));
                
            end
        end
        
        distExpVect(:,1) = parcelDistances(:); 
        distExpVect(:,2) = parcelCoexpression(:); 
        distExpVect(any(isnan(distExpVect), 2), :) = [];
end
%----------------------------------------------------------------------------------
% Fit distance correction according to a defined rule
%----------------------------------------------------------------------------------
%


% select values with distance <100
if strcmp(Fit{1}, 'linear') || strcmp(Fit{1}, 'exp') || strcmp(Fit{1}, 'exp_1_0') || strcmp(Fit{1}, 'decay') || strcmp(Fit{1}, 'exp0') || strcmp(Fit{1}, 'exp1')
    [~,~,c] = GiveMeFit(distExpVect(:,1),distExpVect(:,2),Fit{1});
end

%[param,stat] = sigm_fit(distExpVect(:,1),distExpVect(:,2));

% plot original coexpression-distance .
[xThresholds,yMeans] = BF_PlotQuantiles(distExpVect(:,1),distExpVect(:,2),25,1,1); xlabel('Euclidean distance (mm)'); ylabel('Correlated gene expression');ylim([-1 1]); set(gca,'fontsize',15)
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
    otherwise
        Y = discretize(distExpVect(:,1),xThresholds);
        Residuals = zeros(length(Y),1);
        for val=1:length(Y)
            
            Residuals(val) = distExpVect(val,2) - yMeans(Y(val));
            
        end
        
end

if strcmp(Fit{1}, 'linear') || strcmp(Fit{1}, 'exp') || strcmp(Fit{1}, 'exp_1_0') || strcmp(Fit{1}, 'decay') || strcmp(Fit{1}, 'exp0') || strcmp(Fit{1}, 'exp1')
    hold on; scatter(distExpVect(:,1),FitCurve,1, '.', 'r');
    % get residuals
    Residuals = Rvect - FitCurve;
    BF_PlotQuantiles(distExpVect(:,1),nonzeros(Residuals(:)),50,1,1);   xlabel('Euclidean distance (mm)'); ylabel('Correlated gene expression');ylim([-1 1]); set(gca,'fontsize',15)
else
    BF_PlotQuantiles(distExpVect(:,1),Residuals(:),50,1,1);   xlabel('Euclidean distance (mm)'); ylabel('Correlated gene expression');ylim([-1 1]); set(gca,'fontsize',15)
end


%----------------------------------------------------------------------------------
% Plot corrected sample - sample coexpression matrix
%----------------------------------------------------------------------------------
switch resolution
    case 'sample'
numSamples = size(sampleDistances,1);

% add NaNs to diagonal for reshaping
Idx=linspace(1, size(sampleDistances,1)*size(sampleDistances,1),size(sampleDistances,1));
c=false(1,length(Residuals)+length(Idx));
c(Idx)=true;
nResiduals = nan(size(c));
nResiduals(~c) = Residuals;
nResiduals(c) = NaN;

correctedCoexpression = reshape(nResiduals,[numSamples, numSamples]);
figure; imagesc(correctedCoexpression); caxis([-1,1]);title('Corrected Sample-sample coexpression');
colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);

%----------------------------------------------------------------------------------
% Plot corrected ROI-ROI coexpression matrix
%----------------------------------------------------------------------------------

parcelCoexpression = zeros(length(W),length(W));
parcelDistances = zeros(length(W),length(W));

[sROIs, ind] = sort(ROIs);
correctedCoexpressionSorted = correctedCoexpression(ind, ind);

sampleDistances(logical(eye(size(sampleDistances)))) = NaN;
distancesSorted = sampleDistances(ind, ind);


sampleCoexpression(logical(eye(size(sampleCoexpression)))) = NaN;
coexpressionSorted = sampleCoexpression(ind, ind);

figure; subplot(1,25,[1 2]); imagesc(sROIs);
subplot(1,25,[3 25]); imagesc(coexpressionSorted); caxis([-1,1]); title('Corrected coexpression sorted samples');
colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);

for sub=1:length(W)
    for j=1:length(W)
        
        A = sROIs == W(sub);
        B = sROIs == W(j);
        %for corrected
        if correctDistance == true
            %P = coexpressionSorted(A, B);
            
            P = correctedCoexpressionSorted(A, B);
            
        else
            
            %for uncorrected
            P = coexpressionSorted(A, B);
        end
        D = distancesSorted(A,B);
        parcelDistances(sub,j) = mean(mean(D));
        parcelCoexpression(sub,j) = mean(mean(P));
        
    end
end

end
    
    
    % add zeros values to missing ROIs; p - missing ROI
    
    p = setdiff(nROIs, ROIs);
    expPlot = parcelCoexpression;
    distPlot = parcelDistances;
    if ~isempty(p)
        
        expPlot = insertrows(expPlot,NaN,p);
        expPlot = insertrows(expPlot.', NaN,p).' ; % insert columns
        
        distPlot = insertrows(distPlot,NaN,p);
        distPlot = insertrows(distPlot.', NaN,p).' ; % insert columns
        
    end
    figure; imagesc(expPlot); caxis([-1,1]); title('Parcellation coexpression ROIs');
    colormap([flipud(BF_getcmap('blues',9));[1 1 1];BF_getcmap('reds',9)]);
end
