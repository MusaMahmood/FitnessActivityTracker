function [ Y ] = knnclassify( testSamples, tX, tY, KNN )
%KNNCLASSIFY - Simple K Nearest Neighbor Algorithm using 2norm method:
% Classify using the Nearest neighbor algorithm
% Inputs:
% 	tX          - Train samples
%	tY          - Train labels
%   testSamples - Test  samples
%	KNN		    - Number of nearest neighbors 
%
% Outputs
%	Y           - Predicted targets

    L   = length(tY);
    Uc  = unique(tY);
    if (L < KNN),
       error('You specified more neighbors than there are points.')
    end
    N   = size(testSamples, 1);
    Y= zeros(N,1);
    for i = 1:N,
        dist            = sum( (tX - ones(L,1)*testSamples(i,:) ).^2,2);
        [~, indices]    = sort(dist);  
        n               = hist(tY(indices(1:KNN)), Uc);
        [~, best]       = max(n);
        Y(i)        = Uc(best);
    end
end

