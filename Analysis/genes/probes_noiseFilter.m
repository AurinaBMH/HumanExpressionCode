clear all; close all; 

cd ('data/genes/processedData')
load('MicroarrayDataProbesUpdated.mat')
% select genes that have multiple probes, so thay will be sub-selected for
% comparison
for j=1:2
    % after filtering
    if j==1
        doFilter = true;
    elseif j==2
        doFilter = false;
    end
    
    signalLevel = sum(noiseall,2)./size(noiseall,2);
    indKeepProbes = find(signalLevel>=0.5);
    
    EntrezIDfiltr = DataTableProbe.EntrezID{1,1}(indKeepProbes);
    Expressionallfiltr = Expressionall(indKeepProbes,:);
    
    
    [v, ind] = unique(EntrezIDfiltr);
    
    duplicate_ind = setdiff(1:size(EntrezIDfiltr), ind);
    duplicate_value = unique(EntrezIDfiltr(duplicate_ind));
    
    
    if doFilter
        Ent = EntrezIDfiltr;
        Expr = Expressionallfiltr;
    else
        Ent = DataTableProbe.EntrezID{1,1};
        Expr = Expressionall;
    end
    
    
    for i=1:length(duplicate_value)
        
        
        A = find(Ent==duplicate_value(i));
        %howMany = length(A);
        if length(A)>1
            
            %PL{p} = A;
            r = NaN(length(A));
            for k=1:length(A)
                for l=k+1:length(A)
                    r(k,l) = corr(Expr(A(k),:)', Expr(A(l),:)', 'type', 'Spearman');
                end
                %IND2(w) = A(k);
                
            end
            t=r(:);
            t(isnan(t)) = [];
            
            corVal(i,j) = mean(t);
            
            corVal(i,j) = mean(t);
            
            
            
            %w=w+length(A);
        end
        
    end
end

[p,h,stats] = signrank(corVal(:,1),corVal(:,2));
