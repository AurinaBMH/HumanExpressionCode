%% Author: Aurina
%% Date modified: 2016-03-10
%% Date modified: 2017-07-14
%% This script:
%   1. Loads all microarray data MicroarrayData.mat file
%   2. Uses 1 of 2 possible probe selection options: by max variance or by highest PC
%   3. Finds probes with duplicate entrezIDs
%   4. Chooses one probe from them according to the selected option
%   5. Removes exluded probes from expression data and all other fields.
%   6. Saves MicroarrayDataMaxVAR.mat'or MicroarrayDataPCA.mat with Expressionall, Probe information, Sample information. 
%%

%choose one of those
% To use probes with max variance choose UseMaxVarProbes = 1; keeping the other option = 0;
% To use probes with highest PC choose UseHighestPCProbes = 1; keeping the other option = 0;
% choose if you want to use data with CUST probes or without them:
% UseDataWithCUSTprobes = 1; if UseDataWithCUSTprobes=0, it will load data
% without cust probes. 

UseDataWithCUSTprobes = 0;
ChooseMaxPC = 1;
ChooseMaxVAR = 0;


cd ('/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/');
if UseDataWithCUSTprobes == 1
fprintf(1,'Loading the data with CUST probes and assigning variables\n')
        load('MicroarrayDataWITHCUST.mat');
elseif UseDataWithCUSTprobes == 0
fprintf(1,'Loading the data without CUST probes and assigning variables\n')
        load('MicroarrayData.mat');
end
    

            ProbeName = DataTableProbe.ProbeName{1,1};
            ProbeID = DataTableProbe.ProbeID{1,1};
            EntrezID = DataTableProbe.EntrezID{1,1};
            GeneID = DataTableProbe.GeneID{1,1};
            GeneSymbol = DataTableProbe.GeneSymbol{1,1};
            GeneName = DataTableProbe.GeneName{1,1};
            
                % sigmoid normalise data for MaxVariance calculation        
                ExpressionallNorm = BF_NormalizeMatrix(Expressionall');
                ExpressionallNorm = ExpressionallNorm'; 

            

        
    %% find repeating entrezIDs and calculate variances for them, then take a probe of max variance. 
% 
%         % ------------------------------------------------------------------------------
%         % Find best representative in a set of duplicates using maxVar and remove all others:
%         % ------------------------------------------------------------------------------
        if ChooseMaxVAR == 1
                fprintf(1,'Performing probe selection using Max variance\n')
                Uniq = unique(EntrezID(:,1));
                N = histc(EntrezID, Uniq);
                ProbeList = zeros(length(Uniq),2);

                for k=1:length(Uniq)
                    fprintf(1,'Processing entrez ID %u',k)
                % find indexes for repeating entrexIDs   
                IndForRepEntrezIDs = find(EntrezID==(Uniq(k)));
                    if length(IndForRepEntrezIDs) >=2
                        fprintf(1,'...duplicates found\n')
                  
                
                % take expression values for a selected entrezID
                ExpressionRepEntrezIDs = ExpressionallNorm(IndForRepEntrezIDs,:);
                % calculate variances for expression data for a selected entrezID
                Variance = var(ExpressionRepEntrezIDs,0,2);
                % determine max var value
                [MaxV, IndMaxV] = max(Variance);
                % assign probeID and corresponding entrezID to a list of probes to be used
                % in the analysis
                ProbeList(k,1) = ProbeID(IndForRepEntrezIDs(IndMaxV));
                ProbeList(k,2) = EntrezID(IndForRepEntrezIDs(IndMaxV));
                    else
                ProbeList(k,1) = ProbeID(IndForRepEntrezIDs);
                ProbeList(k,2) = EntrezID(IndForRepEntrezIDs);
                    end

                end
        end
        
         %% PC
        % ------------------------------------------------------------------------------
        % Find best representative in a set of duplicates using PCA and remove all others:
        % ------------------------------------------------------------------------------
        if ChoosePC == 1
        fprintf(1,'Performing probe selection using max PC\n')
                Uniq = unique(EntrezID(:,1));
                N = histc(EntrezID, Uniq);
                ProbeList = zeros(length(Uniq),2);

                for i=1:length(Uniq)
                     fprintf(1,'Processing %u entrez ID\n',i)
                % find indexes for repeating entrexIDs   
                IndForRepEntrezIDs = find(EntrezID==(Uniq(i)));
                    if length(IndForRepEntrezIDs) >=2
                        fprintf(1,'...duplicates found\n')      
                % take expression values for a selected entrezID
                ExpressionRepEntrezIDs = ExpressionallNorm(IndForRepEntrezIDs,:);
                % calculate variances for expression data for a selected
                % entrezID (do not center the data as it is already
                % sigmoided)
                PCAcoeff = pca(ExpressionRepEntrezIDs','Centered',false);
                % determine max PCAcoeff value
                [MaxPCA, IndMaxPCA] = max(PCAcoeff(:,1));
                % assign probeID and corresponding entrezID to a list of probes to be used
                % in the analysis
                ProbeList(i,1) = ProbeID(IndForRepEntrezIDs(IndMaxPCA));
                ProbeList(i,2) = EntrezID(IndForRepEntrezIDs(IndMaxPCA));
                    else
                ProbeList(i,1) = ProbeID(IndForRepEntrezIDs);
                ProbeList(i,2) = EntrezID(IndForRepEntrezIDs);
                    end
                        

                end
        end

                %% check to exclude probes with not maximum variance or max PC and repeating entrezIDs (assign NaN value to
                % them)
                for j=1:length(ProbeID)
                    if ProbeID(j)~=(ProbeList(:,1))
                        ProbeID(j) = NaN;
                        EntrezID(j) = NaN;
                    end
                end

                % combine Probe information variables to a structure
                EntrezID(isnan(EntrezID)) = [];
                ProbeInformation.EntrezID = EntrezID;
                
                GeneID(isnan(ProbeID)) = [];
                ProbeInformation.GeneID = GeneID;
                
                GeneSymbol(isnan(ProbeID)) = [];
                ProbeInformation.GeneSymbol = GeneSymbol;
                
                GeneName(isnan(ProbeID)) = [];
                ProbeInformation.GeneName = GeneName;
                
                ProbeName(isnan(ProbeID)) = [];
                ProbeInformation.ProbeName = ProbeName;
                
                RemoveProbes = ProbeID; 
                
                ProbeID(isnan(ProbeID)) = [];
                ProbeInformation.ProbeID = ProbeID;
                
         for subject=1:6
                % exclude NaN probes keeping 1 probe for 1 entrezID.
                fprintf(1,'Combining and saving the data for subject %u\n', subject)  
                Expression = DataTable.Expression{subject,1};
                Expression(isnan(RemoveProbes),:) = [];
                Expression = Expression';

                
                % combine sample information variables to a structure.
                SampleInformation.StructureNames = DataTable.StructureName{subject,1};
                SampleInformation.MMCoordinates = DataTable.MMcoordinates{subject,1}; 
                SampleInformation.MRIvoxCoordinates = DataTable.MRIvoxCoordinates{subject,1};



                %% save files according to the selected option
                
                if ChoosePC == 1 && UseDataWithCUSTprobes == 0
                
                save(sprintf('MicroarrayDataPCAS0%d.mat', subject), 'Expression', 'ProbeInformation' , 'SampleInformation');
                
                elseif ChooseMaxVAR == 1 && UseDataWithCUSTprobes == 0
                    
                save(sprintf('MicroarrayDataMaxVARS0%d.mat', subject), 'Expression', 'ProbeInformation' , 'SampleInformation');
                
                elseif ChoosePC == 1 && UseDataWithCUSTprobes == 1
                    
                save(sprintf('MicroarrayDataWITHCUSTPCAS0%d.mat', subject), 'Expression', 'ProbeInformation' , 'SampleInformation');
                
                elseif ChooseMaxVAR == 1 && UseDataWithCUSTprobes == 1
                    
                save(sprintf('MicroarrayDataWITHCUSTMaxVARS0%d.mat', subject), 'Expression', 'ProbeInformation' , 'SampleInformation');
                
                end
  %              clearvars -except DataTable DataTableProbe UseMaxVarProbes UseHighestPCProbes UseDataWithCUSTprobes

        end
 
                    
      