function R = compute_reliability(beta,mask,k)
    
    % Compute reliability of betas within a masked region.
    %
    % USAGE: R = compute_reliability(beta,[mask],[k])
    %
    % INPUTS:
    %   beta - [N x 1] cell array, where each cell contains a  matrix of
    %          beta images; K refers to the number of regressors matching
    %          the event name, and V refers to the number of voxels in the mask
    %          image automatically created by SPM
    %   mask (optional) - binary image indicating which voxels to included
    %          in the analysis (default: all voxels)
    %   k (optional) - number of repetitions to train on (default: 1)
    %
    % OUTPUTS:
    %   R - [N x 2] matrix, where R(n,1) is the average within-event
    %       squared error for event n, and R(n,2) is the average between-event squared error
    %
    % Sam Gershman, June 2014
    
    if nargin < 3; mask = 1:size(beta{1},2); end       % use all voxels
    
    mask = mask(mask);
    B = length(beta);
    R = zeros(B,2);
    
    for i = 1:B
        
        y = beta{i}(:,mask);
        N = size(y,1);
        C = nchoosek(1:N,k);
        err = nan(size(C,1),2);
        
        % loop over all combinations
        for j = 1:size(C,1)
            b = mean(y(C(j,:)));    % mean beta estimate
            
            % within-event error
            w = setdiff(1:N,C(j,:));
            if ~isempty(w)
                err(j,1) = mean(mean(bsxfun(@minus,b,y(w,:)).^2));
            end
            
            % between-event error
            for m = 1:B
                if m ~= i
                    err(j,2) = err(j,2) + mean(mean(bsxfun(@minus,b,beta{m}(:,mask)).^2))/(B-1);
                end
            end
        end
        R(i,:) = mean(err);        
    end