function [mask_combined,mask_run1] = fmri_load_combined_mask(EXPT,model,subj)
    
    % Returns a mask obtained by ANDing the masks in all runs
    %
    % USAGE: [mask_combined] = fmri_load_beta(EXPT,model,subj)
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number
    %   subj - subject number
    %
    % OUTPUTS:
    %    
    %   mask_combined - binary image indicating which voxels are included in all run masks
    %   mask_run1 - same for run 1 only, as that is used for complang02_make_langloc_roi
    %
    % Francisco Pereira, September 2014 (adapted from fmri_load_beta)
    % Updates:
    
    S = EXPT.subject(subj);
    M = ['model',num2str(model)];

    mask_combined = [];
    
    for run = 1:length(S.functional)
        % load SPM mask
        rtxt = sprintf('run%d',run);
        V = spm_vol(fullfile(EXPT.analysis_dir,S.name,M,rtxt,'mask.img'));
        mask = spm_read_vols(V); mask = mask~=0;
    
        if run == 1
            mask_combined = mask;
            mask_run1 = mask;
        else
            mask_combined = mask_combined & mask;
        end
    end
