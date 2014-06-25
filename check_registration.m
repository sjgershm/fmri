function check_registration(EXPT,subj,func)
    
    % Check registration of anatomical and functionals.
    %
    % USAGE: check_registration(EXPT,subj)
    % ALT USAGE: check_registration(EXPT,subj,functionalRun)
	
    if nargin < 3, func = 1; end
    S = EXPT.subject(subj);
    niftidir = S.functional(func).niftidir;
    P{1} = fmri_get(fullfile(niftidir,'wf*'));
    P{1} = P{1}(1,:);
    run = S.anatomical.run;
    P{2} = fmri_get(fullfile(S.anatomical.niftidir,sprintf('ws*-%3.4d-*',run)));
    spm_check_registration(char(P));