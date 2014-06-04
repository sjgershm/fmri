function coreg_apply(P,T)
    
    % Apply coregistration transformation to a set of images.
    %
    % USAGE: coreg_apply(P,T)
    %
    % INPUTS:
    %   P - cell array of image names
    %   T - transformation paramgers from spm_coreg
    %
    % Sam Gershman, Jan 2014
    
    M  = spm_matrix(T);
    for j=1:numel(P)
        for k = 1:size(P{j},1)
            MM = spm_get_space(P{j}(k,:));
            spm_get_space(P{j}(k,:),M\MM);
        end
    end