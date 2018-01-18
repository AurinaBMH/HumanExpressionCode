% for each gene calculate average correlation between probes
load('MicroarrayDataProbesUpdated.mat'); 

[v, ind] = unique(DataTableProbe.EntrezID{1});
m=0; w=1; 
for p=1:length(ind)
    A = find(DataTableProbe.EntrezID{1}==v(p)); 
    %howMany = length(A); 
    if length(A)>1
        m=m+1;
        %PL{p} = A;
        r = NaN(length(A));
        for k=1:length(A)
            for l=k+1:length(A)
                r(k,l) = corr(Expressionall(A(k),:)', Expressionall(A(l),:)', 'type', 'Spearman'); 
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

multind = find(corVal); 
corMult = corVal(multind); 

perc = length(find(corMult<0.3))/length(corMult); 
