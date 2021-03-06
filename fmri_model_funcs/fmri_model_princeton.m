function SPM = fmri_model_princeton(EXPT,model,submat,runs)
    
    % First-level GLM analysis.
    %
    % USAGE: SPM = fmri_model(EXPT,model,[submat],[runs])
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number
    %   submat (optional) - vector of subjects to estimate (default: all subjects)
    %   runs (optional) - which runs to analyze
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
                
        S = EXPT.subject(subj); sname = S.name;
        disp(S.name);
        if nargin < 4 || isempty(runs); runs = 1:length(S.functional); end
        
        for i = runs
                        
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
            
            modeldir = fullfile(EXPT.analysis_dir,S.name,['model',num2str(model)],['run',num2str(i)]);
            if exist(modeldir,'dir'); delete(fullfile(modeldir,'*')); end
            mkdir(modeldir);
            cd(modeldir);
            
            % load regressor info (names, onsets and durations)
            para = S.functional(i).para{model};
            if ~isempty(para)
                reg = parse_para(para,EXPT.TR);
                
                % specify image file names
                niftidir = S.functional(i).niftidir;
                run = S.functional(i).run;
                date = S.functional(i).date;
                
                %files{r,1} = fmri_get(fullfile(niftidir,sprintf('f%s_%s-%3.4d-*.nii',sname,date,run)));
                %SPM.xY.P{1} = fmri_get(fullfile(niftidir,sprintf('w*-%3.4d-*',run)));
                SPM.xY.P{1} = fmri_get(fullfile(niftidir,sprintf('wf%s_%s-%3.4d-*.nii',sname,date,run)));
                
                fprintf('run %d\t%d images\n',i,length(SPM.xY.P{1}));
                
                SPM.nscan = size(SPM.xY.P{1},1);
                
                % load movement regressors
                mrp = fullfile(EXPT.analysis_dir,S.name,'movement',['rp',num2str(i)]);
                SPM.Sess.C.C = load(mrp);
                for j = 1:size(SPM.Sess.C.C,2)
                    SPM.Sess.C.name{j} = ['movement',num2str(j)];
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
                        SPM.Sess.U(n) = U;
                    end
                end
                
                % high-pass filter
                SPM.xX.K.HParam = defaults.stats.fmri.hpf;
            end
            
            delete('mask.img'); % make spm re-use directory without prompting
            SPM.xY.P = char(SPM.xY.P);
            SPM = spm_fmri_spm_ui(SPM);
            SPM = spm_spm(SPM);                     %estimate model
            save('SPM','SPM','-v7.3');
            name = SPM.xX.name;
            save('regnames','name');
        end
    end
    
    cd(cdir);       % return to original directory