function [expPlot, parcelCoexpression, correctedCoexpression, Residuals, distExpVect, distPlot] = calculateCoexpression(sampleDistances, selectedGenes, DSvalues, W, ROIs,nROIs, Fit, correctDistance, resolution)

%sampleDistances = maskHalf(sampleDistances);
selectedGenesN = selectedGenes(:,DSvalues(:,1)); % take genes with highest DS values
switch resolution
    case 'sample'
        
        sampleCoexpression = corr(selectedGenesN', 'type', 'Spearman'); % calculate sample-sample coexpression
        %sampleCoexpression = maskHalf(sampleCoexpression);
        
        sampleCoexpression(logical(eye(size(sampleCoexpression)))) = NaN; % replace diagonal with NaN
        distExpVect(:,1) = sampleDistances(:); % make a vector for distances
        distExpVect(:,2) = sampleCoexpression(:); % make a vector for coexpression values
        distExpVect(any(isnan(distExpVect), 2), :) = [];
        
        %distExpVect(isnan(distExpVect,2), : ) = [];  % remove rows for diagonal elememns as thay will have 0 distance and 1 coexpression
        
        Dvect = distExpVect(:,1);
        Rvect = distExpVect(:,2);
        
        figure; imagesc(sampleCoexpression); caxis([-1,1]);title('Sample-sample coexpression');
        colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);
    case 'ROI'
        parcelExpression = zeros(length(W), size(selectedGenesN,2));
        for sub=1:length(W)
            
            A = ROIs == W(sub);
            if length(find(A))>1
                parcelExpression(sub,:) = nanmean(selectedGenesN(A==1,:));
            else
                parcelExpression(sub,:) = selectedGenesN(A==1,:);
            end
            
        end
        
        
        [sROIs, ind] = sort(ROIs);
        distancesSorted = sampleDistances(ind, ind);
        parcelDistances = zeros(length(W),length(W));
        
        for sub=1:length(W)
            for j=1:length(W)
                
                A = sROIs == W(sub);
                B = sROIs == W(j);
                
                D = distancesSorted(A,B);
                parcelDistances(sub,j) = nanmean(D(:));
                
            end
        end
        
        p = setdiff(nROIs, ROIs);
        parcelCoexpression = corr(parcelExpression', 'type', 'Spearman'); % calculate sample-sample coexpression
        parcelCoexpression(logical(eye(size(parcelCoexpression)))) = NaN; % replace diagonal with NaN
        
        if ~isempty(p)
            
            parcelDistances = insertrows(parcelDistances,NaN,p);
            parcelDistances = insertrows(parcelDistances.', NaN,p).' ; % insert columns
            
            parcelCoexpression = insertrows(parcelCoexpression,NaN,p);
            parcelCoexpression = insertrows(parcelCoexpression.', NaN,p).' ; % insert columns
            
        end
        
        distExpVect(:,1) = parcelDistances(:);
        distExpVect(:,2) = parcelCoexpression(:);
        %distExpVect(any(isnan(distExpVect), 2), :) = [];
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
            if ~isnan(distExpVect(val,2))
                Residuals(val) = distExpVect(val,2) - yMeans(Y(val));
            else
                Residuals(val) = NaN;
            end
            
            
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
    case 'ROI'
        numSamples = size(parcelDistances,1);
        sampleDistances = parcelDistances;
        correctedCoexpression = reshape(Residuals,[numSamples, numSamples]);
        
end
figure; imagesc(correctedCoexpression); caxis([-1,1]);title('Corrected coexpression');
colormap([flipud(BF_getcmap('blues',9));[1 1 1];BF_getcmap('reds',9)]);


%----------------------------------------------------------------------------------
% Plot corrected ROI-ROI coexpression matrix
%----------------------------------------------------------------------------------
switch resolution
    case 'sample'
        
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
        
        
        % add NaNs to diagonal for reshaping
        % Idx=linspace(1, size(parcelDistances,1)*size(parcelDistances,1),size(parcelDistances,1));
        % c=false(1,length(Residuals)+length(Idx));
        % c(Idx)=true;
        % nResiduals = nan(size(c));
        % nResiduals(~c) = Residuals;
        % nResiduals(c) = NaN;
        %
        % correctedCoexpression = reshape(nResiduals,[numSamples, numSamples]);
        % figure; imagesc(correctedCoexpression); caxis([-1,1]);title('Corrected Sample-sample coexpression');
        isNormP = zeros(length(W),length(W));
        isNormH = zeros(length(W),length(W));
        
        for sub=1:length(W)
            for j=1:length(W)
                
                A = sROIs == W(sub);
                B = sROIs == W(j);
                %for corrected
                if correctDistance == true
                    %P = coexpressionSorted(A, B);
                    
                    P = correctedCoexpressionSorted(A, B);
%                     if length(P)>4
%                         [isNormH(sub,j),isNormP(sub,j)] = adtest(P(:));
%                     else
%                         isNormH(sub,j) = NaN;
%                         isNormP(sub,j) = NaN;
%                     end
                    
                    
                else
                    
                    %for uncorrected
                    P = coexpressionSorted(A, B);
%                     if length(P)>4
%                         [isNormH(sub,j),isNormP(sub,j)] = adtest(P(:));
%                     else
%                         isNormH(sub,j) = NaN;
%                         isNormP(sub,j) = NaN;
%                     end
                end
                D = distancesSorted(A,B);
                parcelDistances(sub,j) = mean(D(:));
                parcelCoexpression(sub,j) = mean(P(:));
                
            end
        end
        
end

if correctDistance == true && strcmp(resolution, 'ROI')
    expPlot = correctedCoexpression;
else
    expPlot = parcelCoexpression;
end

distPlot = parcelDistances;

if strcmp(resolution, 'sample')
    p = setdiff(nROIs, ROIs);
    if ~isempty(p)
        
        expPlot = insertrows(expPlot,NaN,p);
        expPlot = insertrows(expPlot.', NaN,p).' ; % insert columns
        
        distPlot = insertrows(distPlot,NaN,p);
        distPlot = insertrows(distPlot.', NaN,p).' ; % insert columns
        
    end
end

figure('color','w'); box('off');
imagesc(expPlot); caxis([-1,1]);
if correctDistance
    title('Corrected parcellation coexpression ROIs');
else
    title('Non - corrected parcellation coexpression ROIs');
end
colormap([flipud(BF_getcmap('blues',9));[1 1 1];BF_getcmap('reds',9)]);
end
