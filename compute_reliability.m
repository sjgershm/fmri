function [p, X] = compute_reliability(beta,mask,k,metric)
    
    % Compute reliability of betas within a masked region.
    %
    % USAGE: p = compute_reliability(beta,[mask],[k])
    %
    % INPUTS:
    %   beta - [N x 1] cell array, where each cell contains a  matrix of
    %          beta images; K refers to the number of regressors matching
    %          the event name, and V refers to the number of voxels in the mask
    %          image automatically created by SPM
    %   mask (optional) - binary image indicating which voxels to include
    %          in the analysis (default: all voxels)
    %   k (optional) - number of repetitions to train on (default: 1)
    %
    % OUTPUTS:
    %   R - [N x 2] cell array, where each cell contains a matrix with
    %       nVoxel columns. The number of rows is nchoosek(M,k), where M is the number
    %       repetitions for an event. Each row corresponds to one combination
    %       of k repetitions (the "training set"). R{n,1} is the average within-event
    %       squared error for event in at each voxel (average is taken over repetitions not in the training set);
    %       R{n,2} is the average between-event squared error. If there is
    %       only one repetition, the within-event error will be nan.
    %   p - [N x V] reliability score for each event and voxel. Reliability
    %       is computed as the average difference between between-event and
    %       within-event squared error (thus, positive scores indicate more reliability)
    %
    % Sam Gershman, June 2014
    
    if nargin < 2 || isempty(mask); mask = 1:size(beta{1},2); end       % use all voxels
    if nargin < 3 || isempty(k); k = 1; end
    if nargin < 4 || isempty(metric); metric = 'euclidean'; end
    
    mask = mask(mask);
    B = length(beta);
    p = zeros(B,1);
    X = zeros(B,2);
    
    y = []; n = [];
    for i = 1:B
        y = [y; beta{i}(:,mask)];
        N = size(beta{i},1);
        n = [n; zeros(N,1)+i];
    end
    
    for i = 1:B
        
        ix = find(n==i);
        C = nchoosek(ix,k);
        nc = size(C,1);
        correct = 0;
        count = 0;
        x = zeros(nc,2);
        
        % loop over all combinations
        for j = 1:nc
            b = nanmean(y(C(j,:),:),1);    % mean beta estimate
            w = setdiff(1:length(n),C(j,:));
            nt = n(w);
            r = pdist2(b,y(w,:),metric);
            same = find(nt==i);
            diff = find(nt~=i);
            x(j,:) = [same diff];
            for i1 = 1:length(same)
                for i2 = 1:length(diff)
                    count = count + 1;
                    if r(same(i1)) < r(diff(i2))
                        correct = correct + 1;
                    end
                end
            end
        end
        
        X(i,:) = mean(x);
        p(i) = correct/count;
    end