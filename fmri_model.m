function SPM = fmri_model(EXPT,model,submat)
    
    % First-level GLM analysis.
    %
    % USAGE: SPM = fmri_model(EXPT,model,[submat])
    %
    % INPUTS:
    %   EXPT object
    %   model - model number
    %   submat (optional) - vector of subjects to estimate (default: all subjects)
    %
    % OUTPUTS:
    %   SPM - model structure
    %
    % Sam Gershman, Jan 2014
    
    if nargin < 3 || isempty(submat)
        submat = 1:length(EXPT.subject);
    end
    
    cdir = pwd;
    defaults = spm_get_defaults;
    warning off all
    
    for subj = submat;
        
        % SPM settings
        clear SPM
        SPM.xY.RT = EXPT.TR;
        SPM.xBF.T = defaults.stats.fmri.fmri_t;
        SPM.xBF.T0 = defaults.stats.fmri.fmri_t0;
        SPM.xBF.dt = SPM.xY.RT/SPM.xBF.T;
        SPM.xBF.UNITS   = 'secs';     % time units ('scans', 'secs')
        SPM.xBF.name    = 'hrf';      % basis function type
        SPM.xBF.factor = [];
        SPM.xBF.Volterra = 1;
        SPM.xBF = spm_get_bf(SPM.xBF);
        SPM.xGX.iGXcalc = 'none';     % global intensity normalization (note: 'none' actually means 'session-specific')
        SPM.xVi.form    = 'AR(1)';    % correct for serial correlations ('none', 'AR(1)')
        
        S = EXPT.subject(subj);
        disp(S.name);
        modeldir = fullfile(EXPT.analysis_dir,S.name,['model',num2str(model)]);
        if exist(modeldir,'dir'); delete(fullfile(modeldir,'*')); end
        mkdir(modeldir);
        cd(modeldir);
        
        ii = 0;
        for i = 1:length(S.functional)
            
            % load regressor info (names, onsets and durations)
            para = S.functional(i).para{model};
            if ~isempty(para)
                ii = ii + 1;
                reg = parse_para(para,EXPT.TR);
                
                % specify image file names
                niftidir = S.functional(i).niftidir;
                run = S.functional(i).run;
                SPM.xY.P{ii,1} = fmri_get(fullfile(niftidir,sprintf('w*-%3.4d-*',run)));
                SPM.nscan(ii) = size(SPM.xY.P{ii},1);
                
                % load movement regressors
                mrp = fullfile(EXPT.analysis_dir,S.name,'movement',['rp',num2str(i)]);
                SPM.Sess(ii).C.C = load(mrp);
                for j = 1:size(SPM.Sess(ii).C.C,2)
                    SPM.Sess(ii).C.name{j} = ['movement',num2str(j)];
                end
                
                % configure the input structure array
                n = 0;
                for j=1:numel(reg.onsets)
                    if ~isempty(reg.onsets{j})
                        n = n + 1;
                        U.name = {reg.names{j}};
                        U.ons  = reg.onsets{j}(:);
                        U.dur  = reg.durations(j) .* ones(size(U.ons));
                        U.P    = struct('name', 'none', 'h', 0);
                        SPM.Sess(ii).U(n) = U;
                    end
                end
                
                % high-pass filter
                SPM.xX.K(ii).HParam = defaults.stats.fmri.hpf;
            end
        end
        
        delete('mask.img'); % make spm re-use directory without prompting
        SPM.xY.P = char(SPM.xY.P);
        SPM = spm_fmri_spm_ui(SPM);
        SPM = spm_spm(SPM);                     %estimate model
    end
    
    cd(cdir);       % return to original directory