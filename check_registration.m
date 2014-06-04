function check_registration(EXPT,subj)
    
    % Check registration of anatomical and functionals.
    %
    % USAGE: check_registration(EXPT,subj)
    
    S = EXPT.subject(subj);
    niftidir = S.functional(1).niftidir;
    P{1} = fmri_get(fullfile(niftidir,'wf*'));
    P{1} = P{1}(1,:);
    run = S.anatomical.run;
    P{2} = fmri_get(fullfile(S.anatomical.niftidir,sprintf('ws*-%3.4d-*',run)));
    spm_check_registration(char(P));