function [R,p] = compute_reliability(beta,mask,k)
    
    % Compute reliability of betas within a masked region.
    %
    % USAGE: [R,p] = compute_reliability(beta,[mask],[k])
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
    
    mask = mask(mask);
    V = length(mask);
    B = length(beta);
    R = cell(B,2);
    p = nan(B,V);
    
    for i = 1:B
        
        disp(num2str(i));
        
        y = beta{i}(:,mask);
        N = size(y,1);
        C = nchoosek(1:N,k);
        R{i,1} = zeros(size(C,1),V);
        R{i,2} = R{i,1};
        
        % loop over all combinations
        for j = 1:size(C,1)
            b = nanmean(y(C(j,:)),1);    % mean beta estimate
            
            % within-event error
            w = setdiff(1:N,C(j,:));
            if ~isempty(w)
                R{i,1}(j,:) = nanmean(bsxfun(@minus,b,y(w,:)).^2,1);
            end
            
            % between-event error
            for m = 1:B
                if m ~= i
                    R{i,2}(j,:) = R{i,2}(j,:) + nanmean(bsxfun(@minus,b,beta{m}(:,mask)).^2)/(B-1);
                end
            end
        end
        p(i,:) = nanmean(R{i,2}-R{i,1});
    end