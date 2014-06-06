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
    %          the event name, and V refers to the number of voxels in the mask
    %          image automatically created by SPM
    %   mask - binary image indicating which voxels are included in the mask
    %
    % Sam Gershman, June 2014
    
    S = EXPT.subject(subj);
    M = ['model',num2str(model)];
    load(fullfile(EXPT.analysis_dir,S.name,M,'SPM'));
    
    % load SPM mask
    V = spm_vol(fullfile(EXPT.analysis_dir,S.name,M,'mask.img'));
    mask = spm_read_vols(V); mask = mask~=0;
    
    beta = cell(length(names),1);
    for i = 1:length(names)
        c = find_regressors(SPM.xX.name',names{i});
        c = find(c);
        disp(names{i});
        for j = 1:length(c)
            fname = sprintf('beta_%3.4d.img',c(j));
            V = spm_vol(fullfile(EXPT.analysis_dir,S.name,M,fname));
            Y = spm_read_vols(V);
            beta{i}(j,:) = Y(mask);
        end
    end