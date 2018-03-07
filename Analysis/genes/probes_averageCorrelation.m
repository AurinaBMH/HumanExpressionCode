% for each gene calculate average correlation between probes
cd ('data/genes/processedData')
load('MicroarrayDataWITHcustProbesUpdatedXXX.mat')

signalThreshold = [-1 0.5]; 
figure; colors = [.96 .63 .55; 1 .46 .22]; 
corMult = cell(length(signalThreshold),1);
for i=1:length(signalThreshold)
    
signalLevel = sum(noiseall,2)./size(noiseall,2);
indKeepProbes = find(signalLevel>=signalThreshold(i));


[v, ind] = unique(DataTableProbe.EntrezID{1}(indKeepProbes));
entrezID = DataTableProbe.EntrezID{1}(indKeepProbes); 
Expressionall2 = Expressionall((indKeepProbes),:); 

m=0; w=1; 
corVal = nan(length(ind),1); 
for p=1:length(ind)
    A = find(entrezID==v(p)); 
    %howMany = length(A); 
    if length(A)>1
        m=m+1;
        %PL{p} = A;
        r = NaN(length(A));
        for k=1:length(A)
            for l=k+1:length(A)
                r(k,l) = corr(Expressionall2(A(k),:)', Expressionall2(A(l),:)', 'type', 'Spearman'); 
            end
        %IND2(w) = A(k);
        w=w+1;
        end
        t=r(:); 
        t(isnan(t)) = []; 
        corVal(p) = mean(t);
        %w=w+length(A);
    end
    
end

inds = ~isnan(corVal); 
multind = corVal(inds==1); 
corMult{i} = corVal(multind); 

perc = length(find(corMult{i}<0.3))/length(corMult{i}); 

histogram(corMult{i}, 100,'EdgeColor',[.6 .6 .6],...
    'FaceColor',colors(i,:)); 
xlabel('Average correlation between probes for the same gene','FontSize', 14)
ylabel('Number of genes','FontSize', 14')
set(gcf,'color','w'); hold on; 
end
legendText{1} = sprintf('%d genes (original)', length(corMult{1})); 
legendText{2} = sprintf('%d genes (after QC)', length(corMult{2})); 


legend(legendText{1},legendText{2}); 