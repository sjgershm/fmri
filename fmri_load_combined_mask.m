function [combined_mask] = fmri_load_combined_mask(EXPT,model,subj)
    
    % Returns a mask obtained by ANDing the masks in all runs
    %
    % USAGE: [combined_mask] = fmri_load_beta(EXPT,model,subj)
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number
    %   subj - subject number
    %
    % OUTPUTS:
    %    
    %   combined_mask - binary image indicating which voxels are included in all run masks
    %
    % Francisco Pereira, September 2014 (adapted from fmri_load_beta)
    % Updates:
    
    S = EXPT.subject(subj);
    M = ['model',num2str(model)];

    combined_mask = [];
    
    for run = 1:length(S.functional)
        % load SPM mask
        rtxt = sprintf('run%d',run);
        V = spm_vol(fullfile(EXPT.analysis_dir,S.name,M,rtxt,'mask.img'));
        mask = spm_read_vols(V); mask = mask~=0;
    
        if isempty(combined_mask)
            % this is run 1
            combined_mask = mask;
        else
            combined_mask = combined_mask & mask;
        end
    end
