function masks = fmri_reslice_masks(EXPT,subj,maskdir)
    
    S = EXPT.subject(subj);
    run = S.anatomical.run;
    anatomical = fmri_get(fullfile(S.anatomical.niftidir,sprintf('*-%3.4d-*',run)));
    
    F = dir(fullfile(maskdir,'*.nii'));
    for i = 1:length(F)
        P{i,1} = fullfile(maskdir,F(i).name);
    end
    P = [anatomical; P];
    
    spm_reslice(P);
    
    % load SPM mask
    V = spm_vol(fullfile(EXPT.analysis_dir,S.name,M,'mask.img'));
    mask = spm_read_vols(V); mask = mask~=0;
    
    F = dir(fullfile(maskdir,'r*.nii'));
    for i = 1:length(F)
        f = fullfile(maskdir,F(i).name);
        V = spm_vol(f);
        Y = spm_read_vols(V);
        masks{i} = Y(mask)>0.5;
        delete(f);
    end