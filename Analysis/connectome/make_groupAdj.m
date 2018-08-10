% the script will make
% 1. group matrix
% 2. degree file
% for all the iterations of the pipeline for making plots in neuromarvl.

cd 'data'
load('GenCog_combinations.mat')

ACTs = {'newACT'}; % {'badACT','goodACT', 'noACT'};
seedings = {'dynamic'}; % {'dynamic', 'wm', 'gmwmi'};
siftings = {'SIFT'}; % {'SIFT', 'noSIFT'};
weights = {'DEN'}; % {'NOS', 'FA', 'DEN', 'length'};
parcellations = {'HCP'};
onlyCort = true;
calculateRC = false; 
groupMatrixType = 'consistency';  %{'lengthCV', 'variance', 'consistency'};

% plotting options
scatterColor = [.98 .85 .37];
histoColor = [.28 .75 .57];
edgeColor = [0.45 0.45 0.45];


for i=1:length(ACTs)
    ACT = ACTs{i};
    
    for j = 1:length(seedings)
        %figure;
        %set(gcf,'color','w');
        o=1;
        
        if strcmp(ACT, 'noACT') && j<3
            seeding = seedings{j};
        elseif strcmp(ACT, 'noACT') && j==3
            break
        else
            seeding = seedings{j};
        end
        
        for k = 1:length(siftings)
            sifting = siftings{k};
            
            for l = 1:length(weights)
                weight = weights{l};
                for p=1:length(parcellations)
                    parcel = parcellations{p};
                    
                    matrices = connectomes.(ACT).(seeding).(sifting).(weight).(parcel);
                    matrices = matrices(~cellfun('isempty',matrices));
                    % stack matrices to 3D matrix
                    if onlyCort && strcmp(parcel, 'APARC')
                        LC = 1:34;
                        RC = 42:75;
                        cort = [LC,RC];
                        numNodes = length(cort);
                    elseif onlyCort && strcmp(parcel, 'HCP')
                        LC = 1:180;
                        RC = 191:370;
                        cort = [LC,RC];
                        numNodes = length(cort);
                    elseif onlyCort && strcmp(parcel, 'cust200')
                        LC = 1:100;
                        RC = 111:210;
                        cort = [LC,RC];
                        numNodes = length(cort);
                    else
                        numNodes = size(matrices{1},1);
                    end
                    
                    A = zeros(numNodes, numNodes,length(matrices));
                    for m=1:length(matrices)
                        if onlyCort
                            A(:,:,m) = matrices{m}(cort, cort);
                            matrices{m} = matrices{m}(cort, cort);
                        else
                            A(:,:,m) = matrices{m};
                        end
                    end
                    % make distance matrix: get coordinates
                    coords = connectomes.COG.(parcel)(cort,:);
                    dist = pdist2(coords, coords);
                    % make hemm iid vector
                    hemiid = zeros(numNodes,1);
                    hemiid(1:numNodes/2) = 1;
                    hemiid(numNodes/2+1: numNodes) = 2;
                    
                    % make a group matrix
                    if strcmp(groupMatrixType, 'lengthCV')
                        G = fcn_group_average(A,dist,hemiid);
                    elseif strcmp(groupMatrixType, 'variance')
                        G = giveMeGroupAdj_variance(matrices);
                    elseif strcmp(groupMatrixType, 'consistency')
                        G = giveMeGroupAdj_consistency(matrices, 0.6);
                    end
                    
                    if onlyCort
                        
                        filenameadj = sprintf('adjCort_%s_%s_%s_%s_%s_%s.txt', ACT, seeding, sifting, weight, parcel, groupMatrixType);
                        filenamedeg = sprintf('degCort_%s_%s_%s_%s_%s_%s.txt', ACT, seeding, sifting, weight, parcel, groupMatrixType);
                    else
                        filenameadj = sprintf('adj_%s_%s_%s_%s_%s_%s.txt', ACT, seeding, sifting, weight, parcel, groupMatrixType);
                        filenamedeg = sprintf('deg_%s_%s_%s_%s_%s_%s.txt', ACT, seeding, sifting, weight, parcel, groupMatrixType);
                    end
                    dlmwrite(filenameadj,G,'delimiter','\t');
                    % get the degree distribution
                    deg = degrees_und(G)';
                    
                    %title(sprintf('adj%s%s%s%s%s%s', ACT, seeding, sifting, weight, parcel, groupMatrixType));
                    %o=o+1;
                    
                    dlmwrite(filenamedeg,deg,'delimiter','\t');
                    
                    groupMatr.(ACT).(seeding).(sifting).(weight).(parcel).(groupMatrixType) = G;
                    groupdegr.(ACT).(seeding).(sifting).(weight).(parcel).(groupMatrixType) = deg;
                    
                    % calculate RC curves to determine where the degree
                    % threshold to declare nodes hubs would be
                    
                    
%                     figure;
%                     %subplot(3,2,k);
%                     histogram(deg, 50, 'EdgeColor', edgeColor, 'FaceColor', histoColor);
%                     title (sprintf('degree distribution %s%', thr)); xlabel('Degree')
%                     hold on;
                    
                    % plot the strength distribution of the group matrix
                    %subplot(3,2,k+3);
                    
                    if calculateRC
                        % calculate topological RC
                        numIter = 50;
                        numRepeats = 100;
                        WhatTypeNetwork = 'bu';
                        whatNullModel = 'randmio_und'; %'randmio_und'; %'strength'; %
                        
                        if strcmp(WhatTypeNetwork, 'bu')
                            G = logical(G);
                        end
                        
                        if strcmp(WhatTypeNetwork, 'wuStrength') || strcmp(WhatTypeNetwork, 'wuStrengthBins')
                            kmax = 100;
                        else
                            kmax = max(sum(logical(G)));
                        end
                        strength = sum(G);
                        pThreshold = 0.05;
                        %[~, E] = BF_PlotQuantiles(strength,strength,100);
                        [PhiNorm,PhiTrue,PhiRand] = RichClubPhiNorm(G,kmax, numIter,numRepeats,WhatTypeNetwork,whatNullModel); %, doBins);
                        figure;
                        % Compute p-values
                        pValues = zeros(kmax,1);
                        for w = 1:kmax
                            pValues(w) = mean(PhiTrue(w) <= PhiRand(:,w));
                        end
                        % Significant phi
                        isSig = (pValues <= pThreshold);
                        PhiNormMean = zeros(size(PhiTrue));
                        for w = 1:length(PhiTrue)
                            PhiNormMean(w) = PhiTrue(w)/mean(PhiRand(:,w));
                        end
                        % plot the graphs
                        if strcmp(WhatTypeNetwork, 'wuStrength') || strcmp(WhatTypeNetwork, 'wuStrengthBins')
                            
                            subplot(2,1,1); plot(PhiNormMean, '-','Color','r','LineWidth',2);  xlim([0 100]); ylim([min(PhiNormMean)-0.1 max(PhiNormMean)+0.05]); legend('Phi normalised', 'Location','NorthWest'); ...
                                hold on;
                            plot(find(isSig),PhiNormMean(isSig),'o','Color','r','LineWidth',3); xlim([0 100]); ylim([min(PhiNormMean)-0.1 max(PhiNormMean)+0.05]);
                            title (sprintf('Normalised rich club \n%s', type)); xlabel('Strength'); ylabel('Phi (norm)'); set(gca,'FontSize',12,'fontWeight','bold');
                            subplot(2,1,2); histogram(strength, 50, 'EdgeColor', edgeColor, 'FaceColor', histoColor);  title ('Strength distribution');  xlabel('Strength bins'); set(gca,'FontSize',12,'fontWeight','bold');
                            
                        else
                            
                            subplot(2,1,1); plot(PhiNormMean, '-','Color','r','LineWidth',2);  xlim([min(deg) max(deg)+2]); ylim([0.8 max(PhiNormMean)+0.05]); legend('Phi normalised', 'Location','NorthWest'); ...
                                hold on;
                            plot(find(isSig),PhiNormMean(isSig),'o','Color','r','LineWidth',3); xlim([min(deg) max(deg)+2]); ylim([0.8 max(PhiNormMean)+0.05]);
                            title (sprintf('Normalised rich club')); xlabel('Degree'); ylabel('Phi (norm)'); set(gca,'FontSize',12,'fontWeight','bold');
                            subplot(2,1,2); histogram(deg, 50, 'EdgeColor', edgeColor, 'FaceColor', histoColor); title ('Degree distribution'); xlim([min(deg) max(deg)+2]); xlabel('Degree'); set(gca,'FontSize',12,'fontWeight','bold');
                        end
                    end
                    
                end
            end
        end
    end
end



% make files for bad ACT
% ACTs = {'badACT'};
% seedings = {'dynamic'}; % {'dynamic', 'wm', 'gmwwmi'};
%
%
% for i=1:length(ACTs)
%     ACT = ACTs{i};
%
%     for j = 1:length(seedings)
%
%         seeding = seedings{j};
%     end
%     figure; o=1;
%     set(gcf,'color','w');
%     for k = 1:length(siftings)
%         sifting = siftings{k};
%
%         for l = 1:length(weights)
%             weight = weights{l};
%             for p=1:length(parcellations)
%                 parcel = parcellations{p};
%
%                 matrices = connectomes.(ACT).(seeding).(sifting).(weight).(parcel);
%                 matrices = matrices(~cellfun('isempty',matrices));
%                 % stack matrices to 3D matrix
%                 if onlyCort && strcmp(parcel, 'APARC')
%                     LC = 1:34;
%                     RC = 42:75;
%                     cort = [LC,RC];
%                     numNodes = length(cort);
%                 elseif onlyCort && strcmp(parcel, 'HCP')
%                     LC = 1:180;
%                     RC = 191:370;
%                     cort = [LC,RC];
%                     numNodes = length(cort);
%                 elseif onlyCort && strcmp(parcel, 'cust200')
%                     LC = 1:100;
%                     RC = 111:210;
%                     cort = [LC,RC];
%                     numNodes = length(cort);
%                 else
%                     numNodes = size(matrices{1},1);
%                 end
%
%                 A = zeros(numNodes, numNodes,length(matrices));
%                 for m=1:length(matrices)
%                     if onlyCort
%                         A(:,:,m) = matrices{m}(cort, cort);
%                     else
%                         A(:,:,m) = matrices{m};
%                     end
%                 end
%                 % make distance matrix: get coordinates
%                 coords = connectomes.COG.(parcel)(cort,:);
%                 dist = pdist2(coords, coords);
%                 % make hemm iid vector
%                 hemiid = zeros(numNodes,1);
%                 hemiid(1:numNodes/2) = 1;
%                 hemiid(numNodes/2+1: numNodes) = 2;
%
%                 % make a group matrix
%                 G = fcn_group_average(A,dist,hemiid);
%                 if onlyCort
%
%                     filenameadj = sprintf('adjCort_%s_%s_%s_%s_%s.txt', ACT, seeding, sifting, weight, parcel);
%                     filenamedeg = sprintf('degCort_%s_%s_%s_%s_%s.txt', ACT, seeding, sifting, weight, parcel);
%                 else
%                     filenameadj = sprintf('adj_%s_%s_%s_%s_%s.txt', ACT, seeding, sifting, weight, parcel);
%                     filenamedeg = sprintf('deg_%s_%s_%s_%s_%s.txt', ACT, seeding, sifting, weight, parcel);
%                 end
%                 dlmwrite(filenameadj,G,'delimiter','\t');
%                 % get the degree distribution
%                 deg = degrees_und(G)';
%                 subplot(2,3,o); histogram(deg, 30);
%                 title(sprintf('adj%s%s%s%s%s', ACT, seeding, sifting, weight, parcel));
%                 o=o+1;
%
%                 dlmwrite(filenamedeg,deg,'delimiter','\t');
%
%                 groupMatr.(ACT).(seeding).(sifting).(weight).(parcel) = G;
%                 groupdegr.(ACT).(seeding).(sifting).(weight).(parcel) = deg;
%             end
%         end
%     end
% end





