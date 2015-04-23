function SPM = fmri_model(type,EXPT,model,submat,runs)
    
    % Calls functions that perform first-level GLM analysis.
    %
    % USAGE: SPM = fmri_model(EXPT,model,[submat],[runs])
    %
    % INPUTS:
    %   type: 1-4. default is 1
    %       1: Default HRF modeling, takes 'runs' input for parallelization
    %       2: Princeton version, does take 'runs' input
    %       3: FIR modeling, does not take 'runs' input
    %       4: Old HRF modeling, does not take 'runs' input
    %       5: Single nii HRF modeling, does take 'runs' input
    %   EXPT - experiment structure
    %   model - model number
    %   submat (optional) - vector of subjects to estimate (default: all subjects)
    %   runs (optional) - which runs to analyze
    %
    % OUTPUTS:
    %   SPM - model structure
    %
    % 2015-04-23, Walid Bendris (wbendris@mit.edu): Function created
    
    
    if nargin < 4 || isempty(submat)
        submat = 1:length(EXPT.subject);
    elseif nargin < 5 || isempty(runs)
        runs = [];
    end
    
    switch(type)
        case 1
            SPM = fmri_model_hrf(EXPT,model,submat,runs);
        case 2
            SPM = fmri_model_princeton(EXPT,model,submat,runs);
        case 3
            SPM = fmri_model_fir(EXPT,model,submat);
        case 4
            SPM = fmri_model_old(EXPT,model,submat);
        case 5
            SPM = fmri_model_nii_single(EXPT,model,submat);
    end
            
    
end