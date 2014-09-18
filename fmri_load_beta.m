function [beta, mask] = fmri_load_beta(EXPT,model,subj,names)
    
    % Load beta images.
    %
    % USAGE: [beta, mask] = fmri_load_beta(EXPT,model,subj,names)
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number
    %   subj - subject number
    %   names - [N x 1] cell array of event names
    %
    % OUTPUTS:
    %   beta - [N x 1] cell array, where each cell contains a [K x V] matrix of
    %          beta images; K refers to the number of regressors matching
    %          the event names, and V refers to the number of voxels in the mask
    %          image automatically created by SPM
    %   mask - binary image indicating which voxels are included in the mask
    %
    % Sam Gershman, June 2014
    % Updates:
    %	Walid 2014-08-26: Now can take 'names' as a numeric array to bypass the
    %					  find_regressors step
    
    S = EXPT.subject(subj);
    M = ['model',num2str(model)];
    
    % load SPM mask
    V = spm_vol(fullfile(EXPT.analysis_dir,S.name,M,'run1','mask.img'));
    mask = spm_read_vols(V); mask = mask~=0;
    
    ii = 0;
    for run = 1:length(S.functional)
        load(fullfile(EXPT.analysis_dir,S.name,M,['run',num2str(run)],'regnames'));
        if islogical(names)
            beta = [];
            names = find(names);
            for j = 1:length(names)
                fname = sprintf('beta_%3.4d.img',names(j));
                V = spm_vol(fullfile(EXPT.analysis_dir,S.name,M,fname));
                Y = spm_read_vols(V);
                beta(j,:) = Y(mask);
            end
        else
            beta = cell(length(names),1);
            for i = 1:length(names)
                c = find_regressors(name',names{i});
                c = find(c);
                disp(names{i});
                for j = 1:length(c)
                    fname = sprintf('beta_%3.4d.img',c(j));
                    V = spm_vol(fullfile(EXPT.analysis_dir,S.name,M,fname));
                    Y = spm_read_vols(V);
                    ii = ii + 1;
                    beta{i}(ii,:) = Y(mask);
                end
            end
        end
    end