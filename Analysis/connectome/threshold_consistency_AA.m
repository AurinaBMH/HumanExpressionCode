function [W_thr,Wcv] = threshold_consistency_AA(Ws, p)
%THRESHOLD_CONSISTENCY    Threshold edges ranked by consistency
%
%   W_thr = threshold_consistency(Ws, p);
%
%   This function "thresholds" the group mean connectivity matrix of a set
%   of connectivity matrices by preserving a proportion p (0<p<1) of the
%   edges with the smallest coefficient of variation across the group. All
%   other weights, and all weights on the main diagonal (self-self
%   connections) are set to 0.
%
%   Inputs: Ws,     N-by-N-by-M group of M weighted connectivity matrices
%           p,      proportion of weights to preserve
%                       range:  p=1 (all weights preserved) to
%                               p=0 (no weights removed)
%
%   Output: W_thr,  thresholded group mean connectivity matrix
%
%   Reference: Roberts and Breakspear (2016)
%
%   James Roberts, QIMR Berghofer, 2015
%   Above comments based on BCT's threshold_proportional.m 
%   (BCT: brain-connectivity-toolbox.net)
% 
% can't figure out what to do with the zeros...
% Wlog = log(Ws); % get the log of weights
% Wlog(isinf(Wlog)) = NaN; 
% SDlog = nanstd(Wlog,0,3); 
% SDlog(SDlog==0) = NaN; 

Wcv = zeros(size(SDlog)); 
for i=1:size(SDlog,1)
    for j=i+1:size(SDlog,1)
        Wcv(i,j) = sqrt(exp(SDlog(i,j)^2)-1); % cv for log normal data
    end
end

Wmean=mean(Ws,3); % group mean connectivity matrix
%Wcv=std(Ws,0,3)./Wmean; % coefficient of variation across the group
W_thr=threshold_arbmeasure(Wmean,-Wcv,p); % threshold by low CV = high consistency
