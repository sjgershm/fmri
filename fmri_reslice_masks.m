function masks = fmri_reslice_masks(EXPT,subj,model,maskdir,maskname)
    
    % Reslice NIFTI masks into subject's anatomical space.
    %
    % USAGE: masks = fmri_reslice_masks(EXPT,subj,model,maskdir,[maskname])
    
    S = EXPT.subject(subj);
    run = S.anatomical.run;
    P{1} = fmri_get(fullfile(S.anatomical.niftidir,sprintf('ws*-%3.4d-*',run)));
    
    % load SPM mask
    M = ['model',num2str(model)];
    V = spm_vol(fullfile(EXPT.analysis_dir,S.name,M,'mask.img'));
    mask = spm_read_vols(V); mask = mask~=0;
    
    F = dir(fullfile(maskdir,'*.nii'));
    masks = cell(1,length(F));
    for i = 1:length(F)
        disp(F(i).name);
        P{2} = fullfile(maskdir,F(i).name);
        spm_reslice(P);
        f = fullfile(maskdir,['r',F(i).name]);
        V = spm_vol(f);
        Y = spm_read_vols(V);
        masks{i} = Y(mask)>0.5;
        delete(f);
    end
    
    if nargout > 4
        save(fullfile(EXPT.analysis_dir,S.name,M,maskname),'masks');
    end