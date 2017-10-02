

function [rAdj,numRewirings] = f_shuffleWeights(Adj,numIter);
    % Ben Fulcher, 2014-12-01
    % Shuffles weights, keeping the topology fixed while preserving
    % strenght of each node

    % Get all elements of link data where a connection exists:
%     allActualLinks = Adj(Adj~=0); % Shuffle them 
    [suffled]=random_weight_distribution(Adj);
    
    allActualLinksDataShuffled = suffled;

    % Put them back in the matrix
    rAdj = zeros(size(Adj));
    rAdj(Adj~=0) = allActualLinksDataShuffled;


    % Not relevant to this method:
    numRewirings = 0;
end