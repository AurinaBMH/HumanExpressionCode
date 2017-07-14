 %% combine right/left and cortex/subcortex for each subject separately to make one file for a subject
    subjects = [1 2 3 4 5 6];
    MaxVar = 0;
    PCA = 1;
    Thr = 0;
 Subj = cell(6,1);
    Subj{1} = '/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S01_H0351_2001';
    Subj{2} = '/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S02_H0351_2002';
    Subj{3} = '/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S03_H0351_1009';
    Subj{4} = '/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S04_H0351_1012';
    Subj{5} = '/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S05_H0351_1015';
    Subj{6} = '/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S06_H0351_1016';
%     
%     % choose The parcellation ('aparcaseg', or 'cust100', or 'cust250')
     Parcellations = {'aparcaseg', 'cust100', 'cust250'};
%     % choose the threshold (0, 2, 5, 15 (all points assigned))
     ThresholdList = [2];
load('ProbeInformation.mat');      
   for Parcellation = Parcellations
       for Threshold = ThresholdList
           DataExpression = cell(length(subjects),1);
           DataCoordinates = cell(length(subjects),1);

           
           for subject = subjects
                        S = Subj{subject};
                        if strcmp(Parcellation{1}, 'aparcaseg')
                            NumNodes = 82;
                            Folder = '/default_NativeAnat';
                            FolderName = strcat(S, Folder);
                            cd (FolderName);

                        elseif strcmp(Parcellation{1}, 'cust100')
                            NumNodes = 220;
                            Folder = '/custom100_NativeAnat';
                            FolderName = strcat(S, Folder);
                            cd (FolderName);

                        elseif strcmp(Parcellation{1}, 'cust250')
                            NumNodes = 530;
                            Folder = '/custom250_NativeAnat';
                            FolderName = strcat(S, Folder);
                            cd (FolderName);  
                        end
            
                    
                    
                    if subject == 3 || subject == 4 || subject == 5 || subject == 6
                           IntANDExpressionSorted1 = cell(2,1);
                           Full_information_sorted1 = cell(2,1);
                           load(sprintf('%d_DistThresh%d_S0%d_ExpressionProbePCA_GeneThr%d_leftCortex.mat', NumNodes, Threshold, subject, round(Thr)))
                           IntANDExpressionSorted1{1} = IntANDExpressionSorted;
                           Full_information_sorted1{1} = Full_information_sorted;
                           load(sprintf('%d_DistThresh%d_S0%d_ExpressionProbePCA_GeneThr%d_leftSubcortex.mat', NumNodes, Threshold, subject, round(Thr)))
                           IntANDExpressionSorted1{2} = IntANDExpressionSorted;
                           Full_information_sorted1{2} = Full_information_sorted;
                           
                           IntANDExpression = cat(1, IntANDExpressionSorted1{1},IntANDExpressionSorted1{2});
                           IntANDcoordinates = cat(1, Full_information_sorted1{1},Full_information_sorted1{2});
                    elseif subject == 1 || subject == 2
                           IntANDExpressionSorted1 = cell(4,1);
                           Full_information_sorted1 = cell(4,1);
                           load(sprintf('%d_DistThresh%d_S0%d_ExpressionProbePCA_GeneThr%d_leftCortex.mat', NumNodes, Threshold, subject, round(Thr)))
                           IntANDExpressionSorted1{1} = IntANDExpressionSorted;
                           Full_information_sorted1{1} = Full_information_sorted;
                           load(sprintf('%d_DistThresh%d_S0%d_ExpressionProbePCA_GeneThr%d_leftSubcortex.mat', NumNodes, Threshold, subject,round(Thr) ))
                           IntANDExpressionSorted1{2} = IntANDExpressionSorted;
                           Full_information_sorted1{2} = Full_information_sorted;   
                           load(sprintf('%d_DistThresh%d_S0%d_ExpressionProbePCA_GeneThr%d_rightCortex.mat', NumNodes, Threshold, subject,round(Thr) ))
                           IntANDExpressionSorted1{3} = IntANDExpressionSorted;
                           Full_information_sorted1{3} = Full_information_sorted;
                           load(sprintf('%d_DistThresh%d_S0%d_ExpressionProbePCA_GeneThr%d_rightSubcortex.mat', NumNodes, Threshold, subject,round(Thr) ))
                           IntANDExpressionSorted1{4} = IntANDExpressionSorted;
                           Full_information_sorted1{4} = Full_information_sorted;

                           IntANDExpression = cat(1, IntANDExpressionSorted1{1},IntANDExpressionSorted1{2},IntANDExpressionSorted1{3}, IntANDExpressionSorted1{4});
                           IntANDcoordinates = cat(1, Full_information_sorted1{1},Full_information_sorted1{2},Full_information_sorted1{3}, Full_information_sorted1{4});
                    end
                    SUBJECT = zeros(size(IntANDExpression,1),1);
                    SUBJECT(:,1) = subject;
                   IntANDExpression = sortrows(IntANDExpression,1);
                   IntANDcoordinates = sortrows(IntANDcoordinates,4);
     
                   DataExpression{subject} = [SUBJECT, IntANDExpression];
                   DataCoordinates{subject} = [SUBJECT, IntANDcoordinates];
           end
                   
            
                    
                    %% save the data to S01-S06 combined folder for each parcellation and distance threshold
                 cd ('/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/code/Microarray/S01-S06_combined');
                 fprintf(1, 'Saving combined data for %s parcellation and threshold %u\n', Parcellation{1}, Threshold)
                 if MaxVar == 1
                    save (sprintf('%d_DistThresh%d_%s_combined_ExpressionProbeMaxVar_GeneTher%d.mat', NumNodes, Threshold, Parcellation{1}, round(Thr)), 'DataExpression', 'DataCoordinates', 'ProbeInformation');
                    
                 elseif PCA == 1
                    save (sprintf('%d_DistThresh%d_%s_combined_ExpressionProbePCA_GeneThr%d.mat', NumNodes, Threshold, Parcellation{1}, round(Thr)), 'DataExpression', 'DataCoordinates', 'ProbeInformation');
                    
                 end
      end
   end