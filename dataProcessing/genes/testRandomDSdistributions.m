
% check if distributions are significanly different for genes
pAll = zeros(16202,16202); 
tAll = zeros(16202,16202); 
hAll = zeros(16202,16202); 
for i=1:16202
    for j=i+1:16202
        
    [h,p,ci,stats] = ttest2(DSall(:,i), DSall(:,j)); 
    pAll(i,j) = p; 
    tAll(i,j) = stats.tstat; 
    hAll(i,j) = h; 
    end
    i
end

figure; imagesc(hAll); 

