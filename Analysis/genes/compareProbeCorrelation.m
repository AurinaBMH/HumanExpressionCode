% compare correlation distributions
load('correlations2Probes.mat')
[h,p,ci] = ranksum(r,rall); 
