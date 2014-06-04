function R = compute_reliability(beta,mask)
    
    % Compute reliability of betas within a masked region.
    %
    % USAGE: R = compute_reliability(beta,[mask])
    %
    % INPUTS:
    %   beta - [N x 1] cell array, where each cell contains a  matrix of 
    %          beta images; K refers to the number of regressors matching
    %          the event name, and V refers to the number of voxels in the mask
    %          image automatically created by SPM
    %   mask - binary image indicating which voxels to included in the analysis
    %
    % OUTPUTS:
    %   R - [N x 2] matrix, where R(n,1) is the average within-event
    %   correlation for event n, and R(n,2) is the average correlation
    %   between betas for event n and all other betas.
    %
    % Sam Gershman, June 2014
    
    if nargin < 3; mask = 1:size(beta{1},2); end       % use all voxels
    
    mask = mask(mask);
    B = length(beta);
    Y = []; E = [];
    for i = 1:B
        y = beta{i}(:,mask);
        N = size(y,1);
        Y = [Y; y];
        E = [E; zeros(N,1)+i];
    end
    
    r = corr(Y',Y');    % correlation for every pair of beta images
    r(1:B+1:end)=nan;   % ignore diagonal elements
    
    for i = 1:B
        
        same = E==events(i);
        rs = r(same,same);      % within-event correlation
        diff = E~=events(i);
        rd = r(diff,diff);      % between-event correlation
        
        R(i,:) = [nanmean(rs(:)) nanmean(rd(:))];
        
    end