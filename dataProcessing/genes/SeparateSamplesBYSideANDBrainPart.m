%% Author: Aurina
%Last modiffied: 2016-04-15

%This script loads MicroarrayDataPCAS[subject].mat file and separates samples according to Allen labels to right/left, cortex/subcortex
%resulting in 4 files for each subject. 
%% choose the options: subjects (subjects [ 1 2 3 4 5 6]), side (KeepRight = 1, KeepLeft = 0), brain part (KeepCortex = 1/0)


clear all; close all; 
% for subjects 1 and 2 do both right and left' for subjects 3:6 do only
% left (run script for those separately);

Thr = 0;
KeepCortex = 0 ; % choose =1 (if you want to keep cortical samples; choose =0 if you want to keep subcortical samples)
KeepRight = 0; % 1 - keep right, 0 -not keep rightchoose to keep right or left (one or the other)

if KeepRight == 1
    subjects = [1 2];
elseif KeepRight == 0
    subjects = [1 2 3 4 5 6];
end


%% Separate cortex/subcortex and right/left samples according to Allen institute Structure names. 

for subject = subjects
    % load data file
     FileName = sprintf('MicroarrayDataPCAS0%d_GeneThr%d.mat', subject, Thr);

    load(FileName);
    Cortex = {'gyrus', 'gyri', 'sulcus', 'occipital pole', 'planum temporale', 'cuneus', ... 
        'frontal pole', 'operculum', 'planum polare', 'temporal pole' , 'paracentral lobule', 'cortex'};
    
    % Keep only cortex samples
    index = cell(length(Cortex),1);

    for i=1:length(SampleInformation.StructureNames)
        Structure = SampleInformation.StructureNames{i};
        %check if structure name has a part related to cortex
        for j = 1:length(Cortex)
            index{j} = findstr(Structure, Cortex{j});
        end
        %check if cortex related name part was found
        emptyCells = cellfun(@isempty,index);
    


            if KeepCortex == 1 % will exclude subcortical samples
                
                    BrainPart = {'Cortex'};
                    if sum(emptyCells) == length(Cortex)   
                        Expression(i,:) = NaN;
                        SampleInformation.MRIvoxCoordinates(i,:) = NaN;
                        SampleInformation.MMCoordinates(i,:) = NaN;
                        SampleInformation.StructureNames{i} = NaN;
                    end
                    
            elseif KeepCortex == 0  % will exclude cortical samples
                
                    BrainPart = {'Subcortex'};
                    if sum(emptyCells) < length(Cortex)
                        Expression(i,:) = NaN;
                        SampleInformation.MRIvoxCoordinates(i,:) = NaN;
                        SampleInformation.MMCoordinates(i,:) = NaN;
                        SampleInformation.StructureNames{i} = NaN;

                    end
            end
    end

        
 

    % exclude nonexisting data
     Expression(any(isnan(Expression),2),:) = [];
     SampleInformation.MRIvoxCoordinates(any(isnan(SampleInformation.MRIvoxCoordinates),2),:) = [];
     SampleInformation.MMCoordinates(any(isnan(SampleInformation.MMCoordinates),2),:) = [];
     SampleInformation.StructureNames=SampleInformation.StructureNames(cell2mat(cellfun(@ischar,SampleInformation.StructureNames,'UniformOutput',0)));
     
     %% separate samples by side: left/right
     Sides = {'left', 'right'};

      for k=1:length(SampleInformation.StructureNames)
        Structure2 = SampleInformation.StructureNames{k};
        if KeepRight == 1
            Side = Sides{2};
            index2 = findstr(Structure2, Side);
        elseif KeepRight == 0
            Side = Sides{1};
            index2 = findstr(Structure2, Side);
        end
        % exclude not relevant side samples
        if isempty(index2)
            Expression(k,:) = NaN;
            SampleInformation.MRIvoxCoordinates(k,:) = NaN;
            SampleInformation.MMCoordinates(k,:) = NaN;
            SampleInformation.StructureNames{k} = NaN;
        end
            
      end
     
     % delete NaN rows.
     Expression(any(isnan(Expression),2),:) = [];
     SampleInformation.MRIvoxCoordinates(any(isnan(SampleInformation.MRIvoxCoordinates),2),:) = [];
     SampleInformation.MMCoordinates(any(isnan(SampleInformation.MMCoordinates),2),:) = [];
     SampleInformation.StructureNames=SampleInformation.StructureNames(cell2mat(cellfun(@ischar,SampleInformation.StructureNames,'UniformOutput',0)));
     %% save data to file


    
     SaveFileName = sprintf('MicroarrayDataPCAS0%d_GeneThr%d_%s%s.mat', subject, Thr, Side, BrainPart{1});
     save (SaveFileName, 'Expression', 'ProbeInformation', 'SampleInformation');

    end



%%
% 
% subjects = [1 2 3 4 5 6];
% 
% KeepCortex = 1; % choose =1 (if you want to keep cortical samples; choose =0 if you want to keep subcortical samples)
% KeepRight = [0 1]; % choose to keep right or left (one or the other)
% KeepLeft = 0; 
% 
% Thr = 3;
% %% Separate cortex/subcortex and right/left samples according to Allen institute Structure names. 
% 
% for subject = subjects
%     % load data file
%      FileName = sprintf('MicroarrayDataPCAS0%d_GeneThr%d.mat', subject, Thr);
% 
%     load(FileName);
%     Cortex = {'gyrus', 'gyri', 'sulcus', 'occipital pole', 'planum temporale', 'cuneus', ... 
%         'frontal pole', 'operculum', 'planum polare', 'temporal pole' , 'paracentral lobule', 'cortex'};
%     
%     % Keep only cortex samples
%     index = cell(length(Cortex),1);
% 
%     for i=1:length(SampleInformation.StructureNames)
%         Structure = SampleInformation.StructureNames{i};
%         %check if structure name has a part related to cortex
%         for j = 1:length(Cortex)
%             index{j} = findstr(Structure, Cortex{j});
%         end
%         %check if cortex related name part was found
%         emptyCells = cellfun(@isempty,index);
%         
%             if KeepCortex == 1 % will exclude subcortical samples
%                 
%                     BrainPart = {'Cortex'};
%                     if sum(emptyCells) == length(Cortex)   
%                         Expression(i,:) = NaN;
%                         SampleInformation.MRIvoxCoordinates(i,:) = NaN;
%                         SampleInformation.MMCoordinates(i,:) = NaN;
%                         SampleInformation.StructureNames{i} = NaN;
%                     end
%             elseif KeepCortex == 0 % will exclude cortical samples
%                 
%                     BrainPart = {'Subcortex'};
%                     if sum(emptyCells) < length(Cortex)
%                         Expression(i,:) = NaN;
%                         SampleInformation.MRIvoxCoordinates(i,:) = NaN;
%                         SampleInformation.MMCoordinates(i,:) = NaN;
%                         SampleInformation.StructureNames{i} = NaN;
% 
%                     end
%             end
% 
%     end 
%  
% 
%     % exclude nonexisting data
%      Expression(any(isnan(Expression),2),:) = [];
%      SampleInformation.MRIvoxCoordinates(any(isnan(SampleInformation.MRIvoxCoordinates),2),:) = [];
%      SampleInformation.MMCoordinates(any(isnan(SampleInformation.MMCoordinates),2),:) = [];
%      SampleInformation.StructureNames=SampleInformation.StructureNames(cell2mat(cellfun(@ischar,SampleInformation.StructureNames,'UniformOutput',0)));
%      
%      %% separate samples by side: left/right
%      Sides = {'left', 'right'};
% 
%       for k=1:length(SampleInformation.StructureNames)
%         Structure2 = SampleInformation.StructureNames{k};
%         if KeepRight == 1
%             Side = Sides{2};
%             index2 = findstr(Structure2, Side);
%         elseif KeepRight == 0
%             Side = Sides{1};
%             index2 = findstr(Structure2, Side);
%         end
%         % exclude not relevant side samples
%         if isempty(index2)
%             Expression(k,:) = NaN;
%             SampleInformation.MRIvoxCoordinates(k,:) = NaN;
%             SampleInformation.MMCoordinates(k,:) = NaN;
%             SampleInformation.StructureNames{k} = NaN;
%         end
%             
%       end
%      
%      % delete NaN rows.
%      Expression(any(isnan(Expression),2),:) = [];
%      SampleInformation.MRIvoxCoordinates(any(isnan(SampleInformation.MRIvoxCoordinates),2),:) = [];
%      SampleInformation.MMCoordinates(any(isnan(SampleInformation.MMCoordinates),2),:) = [];
%      SampleInformation.StructureNames=SampleInformation.StructureNames(cell2mat(cellfun(@ischar,SampleInformation.StructureNames,'UniformOutput',0)));
%      %% save data to file
%      
%      SaveFileName = sprintf('MicroarrayDataPCAS0%d_GeneThr%d_%s%s.mat', subject, Thr, Side, BrainPart{1});
%      save (SaveFileName, 'Expression', 'ProbeInformation', 'SampleInformation');
%     
% end

