function fmri_preproc_princeton(EXPT,subj,tasks)
    
    % Preprocess fMRI data.
    %
    % USAGE: fmri_preproc(EXPT,subj,[tasks])
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   subj - subject number
    %   tasks (optional) - if one of the following strings, it will perform
    %                      the specified preprocessing step:
    %                       'dicom_convert' - convert dicoms to nifti
    %                       'realign' - motion correction
    %                       'coregister' - image coregistration
    %                       'normalize' - warp anatomical to template
    %                       'smooth' - smooth images with gaussian kernel
    %                       'all' - all of the above (in order) except smoothing
    %                      if tasks is a cell array, it will do the
    %                      specified tasks in sequence. By default, tasks = 'all'
    %
    % OUTPUTS:
    %   Dicom conversion writes nifti files to 'nii' directory.
    %   Realignment and coregistration induce changes in the headers of the
    %   nifti files. Normalization writes out new files with 'w' prefix.
    %   Smoothing writes out new files with 's' prefix.
    %
    % Sam Gershman, Jan 2014
    
    
    if nargin < 3 || isempty(tasks); tasks = 'all'; end
    
    if iscell(tasks)
        for i = 1:length(tasks)
            fmri_preproc(EXPT,subj,tasks{i});
        end
        return
    end
    
    if isfield(EXPT,'multiple_anatomicals')
        multiple_anatomicals = EXPT.multiple_anatomicals;
    else
        multiple_anatomicals = 0;
    end
    sessions = unique([EXPT.subject(1).functional(:).session]);
    nSessions = length(sessions);
        
    S = EXPT.subject(subj);
    adir = fullfile(EXPT.analysis_dir,S.name);  % make analysis directory
    if ~exist(adir,'dir'); mkdir(adir); end
    
    switch tasks
      case 'all'
        % Do everything
        tasks = {'dicom_convert' 'realign' 'coregister' 'normalize'};
        fmri_preproc(EXPT,subj,tasks);
        return
        
      case 'dicom_convert'
        % Convert dicom files to nifti.
        
        disp('Converting dicoms to nifti...');
        
        curdir = pwd;            

        for r = 1:length(S.functional)
            %disp(['run ',num2str(r)]);
            dicomdir = S.functional(r).dicomdir;
            niftidir = S.functional(r).niftidir;
            run = S.functional(r).run;
            % FP:
            %files = dir(fullfile(dicomdir,sprintf('*-%d-*',run)));
            files = dir(fullfile(dicomdir,sprintf('%d-*.dcm',run)));
            fprintf('run %d\t%d files\n',r,length(files));

            files = dir2char(files,dicomdir);
            hdr = spm_dicom_headers(files);
            if ~exist(niftidir,'dir'); mkdir(niftidir);  end
            delete(fullfile(niftidir,'w*'));
            delete(fullfile(niftidir,'s*'));
            cd(niftidir);
            spm_dicom_convert(hdr,'all','flat','nii');
        end

        disp('anatomical');
        if multiple_anatomicals
            for ia = 1:length(S.anatomical)
                dicomdir = S.anatomical(ia).dicomdir;
                niftidir = S.anatomical(ia).niftidir;
                run = S.anatomical(ia).run;
                %FP
                %files = dir(fullfile(dicomdir,sprintf('*-%d-*',run)));
                files = dir(fullfile(dicomdir,sprintf('%d-*.dcm',run)));
                files = dir2char(files,dicomdir);
            
                hdr = spm_dicom_headers(files);
                if ~exist(niftidir,'dir'); mkdir(niftidir); end
                cd(niftidir);
                spm_dicom_convert(hdr,'all','flat','nii');
                cd(curdir);
            end
        else
            dicomdir = S.anatomical.dicomdir;
            niftidir = S.anatomical.niftidir;
            run = S.anatomical.run;
            %FP
            %files = dir(fullfile(dicomdir,sprintf('*-%d-*',run)));
            files = dir(fullfile(dicomdir,sprintf('%d-*.dcm',run)));
            files = dir2char(files,dicomdir);
            
            hdr = spm_dicom_headers(files);
            if ~exist(niftidir,'dir'); mkdir(niftidir); end
            cd(niftidir);
            spm_dicom_convert(hdr,'all','flat','nii');
            cd(curdir);
        end
        
      case 'realign'
        % This step does motion correction. Only a mean image (mean_*) is
        % resliced; all other images have their headers modified.
        
        disp('Realigning...')
        
        if multiple_anatomicals == 0
            
            % get nifti filenames
            sname = EXPT.subject(subj).name;
            
            for r = 1:length(S.functional)
                niftidir = S.functional(r).niftidir;
                run = S.functional(r).run;
                date = S.functional(r).date;
                
                files{r,1} = fmri_get(fullfile(niftidir,sprintf('f%s_%s-%3.4d-*.nii',sname,date,run)));
                fprintf('run %d\t%d files\n',r,length(files{r,1}));
                %for e = 1:length(files{r,1})
                %    fprintf('\t%s\n',files{r,1}(e,:));
                %end
            end
            
            spm_realign(files); % run realignment
            spm_reslice(files,struct('mean',1,'which',0));  % write mean image
            
            % move the movement parameter files (rp*) to analysis directory
            E = fullfile(EXPT.analysis_dir,S.name,'movement');
            if ~exist(E,'dir'); mkdir(E); end
            for r = 1:length(S.functional)
                niftidir = S.functional(r).niftidir;
                run = S.functional(r).run;
                rp = fmri_get(fullfile(niftidir,sprintf('rp*-%3.4d-*',run)));
                movefile(rp,fullfile(E,['rp',num2str(r)]));
            end

        else; % multiple anatomicals -------------------------------
            
            % get nifti filenames
            sname = EXPT.subject(subj).name;
            
            for s = sessions
                fprintf('\tsession %d\n',s);
            
                files = {};
                ridx = 1;
                for r = 1:length(S.functional)
                    niftidir = S.functional(r).niftidir;
                    run = S.functional(r).run;
                    date = S.functional(r).date;
                    session = S.functional(r).session;
                    
                    if s == session
                        files{ridx,1} = fmri_get(fullfile(niftidir,sprintf('f%s_%s-%3.4d-*.nii',sname,date,run)));
                        fprintf('\t\trun %d\t%d files\n',r,length(files{ridx,1}));
                        ridx = ridx + 1;
                    end
                    %for e = 1:length(files{r,1})
                    %    fprintf('\t%s\n',files{r,1}(e,:));
                    %end
                end
            
                spm_realign(files); % run realignment
                spm_reslice(files,struct('mean',1,'which',0));  % write mean image
            
                % move the movement parameter files (rp*) to analysis directory
                %E = fullfile(EXPT.analysis_dir,S.name,'movement');
                E = fullfile(EXPT.analysis_dir,S.name,S.dir_movement);
                
                if ~exist(E,'dir'); mkdir(E); end
                for r = 1:length(S.functional)
                    niftidir = S.functional(r).niftidir;
                    run = S.functional(r).run;
                    date = S.functional(r).date;
                    session = S.functional(r).session;
                    
                    if s == session
                        rp = fmri_get(fullfile(niftidir,sprintf('rp*-%3.4d-*',run)));
                        movefile(rp,fullfile(E,['rp',num2str(r)]));
                    end
                end
            end; % for over sessions
        end
  
        
      case 'coregister'
        % This step first coregisters the mean functional to the
        % anatomical, and then coregisters the anatomical to the
        % MNI template.
        
        disp('Coregistration...')
        sname = EXPT.subject(subj).name;

        if multiple_anatomicals == 0
        
            % 1) mean functional -> anatomical
            % use normalied mutual information for registering images from
            % different modalities (epi -> T1)
            niftidir = S.functional(1).niftidir;
            mean_epi = fmri_get(fullfile(niftidir,'mean*'))
            run = S.anatomical.run;
            %anatomical = fmri_get(fullfile(S.anatomical.niftidir,sprintf('*-%3.4d-*',run)))
            anatomical = fmri_get(fullfile(S.anatomical.niftidir,sprintf('*-%3.4d-*',run)))
            %pause
            
            T1 = spm_coreg(anatomical,mean_epi,struct('cost_fun','nmi','graphics',0));
            
            % transform all other functionals

            for r = 1:length(S.functional)
                niftidir = S.functional(r).niftidir;
                run = S.functional(r).run;
                date = S.functional(r).date;
                session = S.functional(r).session;
                
                %P{r,1} = fmri_get(fullfile(niftidir,sprintf('f*-%3.4d-*',run)));
                P{r,1} = fmri_get(fullfile(niftidir,sprintf('f%s_%s-%3.4d-*.nii',sname,date,run)));
                
                fprintf('run %d\t%d files\n',r,length(P{r,1}));
                %for e = 1:length(P{r,1})
                %    fprintf('\t%s\n',P{r,1}(e,:));
                %end
            end
            %pause
            
            P = [P; mean_epi];
            coreg_apply(P,T1);
            
            % 2) anatomical -> MNI template
            % use normalied cross-correlation for registering images from
            % the same modality (T1 -> T1)
            template = which('T1.nii');
            T2 = spm_coreg(template,anatomical,struct('cost_fun','ncc','graphics',0));
            P = [P; anatomical];
            coreg_apply(P,T2);   % transform all other functionals and anatomical
            
            % save coregistration parameters
            save(fullfile(adir,'coreg_params'),'T1','T2');
        
        else; % multiple anatomicals -------------------------------     
            
            for s = sessions
                fprintf('\tsession %d\n',s);
                
                % 1) mean functional -> anatomical
                % use normalied mutual information for registering images from
                % different modalities (epi -> T1)
                niftidir = S.anatomical(s).niftidir;
                mean_epi = fmri_get(fullfile(niftidir,'mean*'))
                run = S.anatomical(s).run;
                anatomical = fmri_get(fullfile(S.anatomical(s).niftidir,sprintf('*-%3.4d-*',run)))
                %pause
                
                T1 = spm_coreg(anatomical,mean_epi,struct('cost_fun','nmi','graphics',0));
                
                ridx = 1;
                P = {};
                
                % transform all other functionals
                for r = 1:length(S.functional)
                    niftidir = S.functional(r).niftidir;
                    run = S.functional(r).run;
                    date = S.functional(r).date;
                    session = S.functional(r).session;
                    
                    if s == session
                        %P{r,1} = fmri_get(fullfile(niftidir,sprintf('f*-%3.4d-*',run)));
                        P{ridx,1} = fmri_get(fullfile(niftidir,sprintf('f%s_%s-%3.4d-*.nii',sname,date,run)));
                        fprintf('\t\trun %d\t%d files\n',r,length(P{ridx,1}));
                        %for e = 1:length(P{r,1})
                        %    fprintf('\t%s\n',P{r,1}(e,:));
                        %end
                        ridx = ridx + 1;
                    end
                end
                %pause
                
                P = [P; mean_epi];
                coreg_apply(P,T1);
                
                % 2) anatomical -> MNI template
                % use normalied cross-correlation for registering images from
                % the same modality (T1 -> T1)
                template = which('T1.nii');
                T2 = spm_coreg(template,anatomical,struct('cost_fun','ncc','graphics',0));
                P = [P; anatomical];
                coreg_apply(P,T2);   % transform all other functionals and anatomical
                
                % save coregistration parameters
                save(fullfile(adir,sprintf('coreg_params_session%d',s)),'T1','T2'); 
                
            end; % of loop over sessions
        end; % of if on multiple anatomicals
 
        
      case 'normalize'
        % Segment and normalize anatomical image to MNI template.
        % Functionals and anatomical prepended with 'w'.
        
        disp('Normalizing anatomical...')
        sname = EXPT.subject(subj).name;
        
        if multiple_anatomicals == 0
 
            run = S.anatomical.run
            anatomical = fmri_get(fullfile(S.anatomical.niftidir,sprintf('*-%3.4d-*',run)))
            
            res = spm_preproc(anatomical);   % compute warping parameters
            sn = spm_prep2sn(res);
            spm_write_sn(anatomical,sn);     % normalize anatomical using the warp parameters already calculated
            
            % normalize functionals
            for r = 1:length(S.functional)
                niftidir = S.functional(r).niftidir;
                run = S.functional(r).run;
                date = S.functional(r).date;
                
                %P{r,1} = fmri_get(fullfile(niftidir,sprintf('f%s_%s-%3.4d-*.nii',sname,date,run)));
                
                %P = fmri_get(fullfile(niftidir,sprintf('f*-%3.4d-*',run)));
                P = fmri_get(fullfile(niftidir,sprintf('f%s_%s-%3.4d-*.nii',sname,date,run)));
                
                fprintf('run %d\t%d files\n',r,length(P));
                %for e = 1:length(P)
                %    fprintf('\t%s\n',P(e,:));
                %end
                %pause
                
                spm_write_sn(P,sn);
            end
                        
            adir = fullfile(EXPT.analysis_dir,EXPT.subject(subj).name);
            if ~exist(adir,'dir'); mkdir(adir); end
            save(fullfile(adir,'normalization_params'),'sn'); 
                
        else; % multiple anatomicals -------------------------------     
            
            for s = sessions
                fprintf('\tsession %d\n',s);
            
                run = S.anatomical(s).run
                anatomical = fmri_get(fullfile(S.anatomical(s).niftidir,sprintf('s*-%3.4d-*',run)))
                niftidir = S.anatomical(s).niftidir;
                mean_epi = fmri_get(fullfile(niftidir,'mean*'))
                
                res = spm_preproc(anatomical);   % compute warping parameters
                sn = spm_prep2sn(res);
                %%spm_write_sn(anatomical,sn);     % normalize anatomical using the warp parameters already calculated
            
                % normalize functionals
                for r = 1:length(S.functional)
                    niftidir = S.functional(r).niftidir;
                    run = S.functional(r).run;
                    date = S.functional(r).date;
                    session = S.functional(r).session;

                    if s == session
                        P = fmri_get(fullfile(niftidir,sprintf('f%s_%s-%3.4d-*.nii',sname,date,run)));
                        fprintf('run %d\t%d files\n',r,length(P));
                        %%spm_write_sn(P,sn);
                    end
                end

                % also do it to the mean epi
                spm_write_sn(mean_epi,sn);
 
                %%adir = fullfile(EXPT.analysis_dir,EXPT.subject(subj).name);
                %%if ~exist(adir,'dir'); mkdir(adir); end
                %%save(fullfile(adir,sprintf('normalization_params_session%d',s)),'sn');
                            
            end; % of for loop over sessions
        end
        
            
      case 'smooth'
        % Smooth images with Gaussian kernel (width specified in EXPT.fwhm).
        % Functionals prepended with 's'.
        
        disp('Smoothing...');
        
        for r = 1:length(S.functional)
            niftidir = S.functional(r).niftidir;
            run = S.functional(r).run;
            P = fmri_get(fullfile(niftidir,sprintf('w*-%3.4d-*',run)));
            for j = 1:size(P,1)
                [pth,nam,ext,num] = spm_fileparts(P(j,:));
                u = fullfile(pth,['s' nam ext num]);    % preprend 's' to filenames
                spm_smooth(P(j,:),u,EXPT.fwhm);
            end
        end


        
      case 'createDatasetMPSC'
        % create a MATLAB dataset with post-processed data and per-TR labels
        % (run AFTER fmri_model_princeton)

        model = 3;
        M = ['model',num2str(model)];
        
        fprintf('extracting event information\n');
        
        %% generate regressors per run
        
        runsToUse = ones(length(S.functional),1);
        
        ridx = 1;
        for s = sessions
            fprintf('\tsession %d\n',s);
            
            for r = 1:length(S.functional)
                niftidir = S.functional(r).niftidir;
                rname    = S.functional(r).name;
                run      = S.functional(r).run;
                date     = S.functional(r).date;
                session  = S.functional(r).session;
                rlength  = S.functional(r).length; % in TRs
                
                switch rname
                  case {'stim_sentencesAllSubset_1',...
                        'stim_sentencesAllSubset_2',...
                        'stim_sentencesAllSubset_3',...
                        'stim_sentencesAllSubset_4',...
                        'stim_sentencesAllSubset_5',...
                        'stim_sentencesAllSubset_6',...
                        'stim_sentencesAllSubset_7',...
                        'stim_sentencesAllSubset_8'}
                    runsToUse(r) = 1;
                  case {'stim_sentencesSubsetA_1';
                        'stim_sentencesSubsetB_1'}
                    runsToUse(r) = 1;
                  case {'stim_words1',...
                        'stim_words2',...
                        'stim_words3',...
                        'stim_words4',...
                        'stim_words5',...
                        'stim_words6',...
                        'stim_words7',...
                        'stim_words8'}
                    runsToUse(r) = 1;
                  otherwise
                    runsToUse(r) = 0;
                end
                    
                if (s == session) & (runsToUse(r) == 1)
                    % want the one with event durations
                    parafile = S.functional(r).para{1};
                    
                    % get all the events/durations in units of TRs 
                    P = parse_para(parafile);
                    nEventTypes = length(P.events); % in this run
                    nEventMax   = length(P.names_all);
                    
                    % event regressors
                    %eventRegressors{r}
                    tmp = zeros(rlength,nEventMax);

                    fprintf('\t\trun %d - %s\n',r,rname);
                    for e = 1:nEventTypes
                        index    = P.events(e); % in the names_all list
                        duration = P.durations(e); % always the same
                        onsets   = P.onsets{e}; % in this run
                        
                        %fprintf('\t\t\t%s\t%1.2f\t%d\n',P.names_all{index},duration,length(onsets));
                        
                        for o = 1:length(onsets)
                            onset = onsets(o);
                            
                            if round(onset) == onset
                                % precisely on TR boundary
                                TR = onset + 1;
                                durationInTR = 1;
                            else
                                TR = ceil(onset);
                                durationInTR = TR-onset;
                            end
                            tmp(TR,index) = durationInTR;
                            
                            leftover = duration - durationInTR;
                            while leftover > 0
                                TR = TR + 1;
                                if leftover >= 1
                                    tmp(TR,index) = 1;
                                else
                                    tmp(TR,index) = leftover;
                                end
                                leftover = leftover - 1;
                            end
                        end                    
                    end; % for loop over events

                    eventRegressorsPerRun{r} = tmp;
                    otherRegressorsPerRun{r} = [repmat(ridx,rlength,1),repmat(s,rlength,1)];
                    
                    %imagesc(tmp(:,P.events)); pause
                                        
                    %P = fmri_get(fullfile(niftidir,sprintf('f%s_%s-%3.4d-*.nii',sname,date,run)));
                    %fprintf('run %d\t%d files\n',r,length(P));
                    
                    ridx = ridx + 1; % index over runs actually used
                end; % if for runs in this session
            end; % loop over runs
        end; % loop over sessions
        
        %% aggregate all the regressors across runs
        
        ntotal = 0;
        nmax = size(eventRegressorsPerRun{1},2);
        for r = 1:length(S.functional)
            if runsToUse(r)
                if size(eventRegressorsPerRun{r},2) ~= nmax
                    fprintf('error: mismatch in #events\n');return;
                end
                ntotal = ntotal + size(eventRegressorsPerRun{r},1);
            end
        end
            
        eventRegressors = zeros(ntotal,nmax);
        otherRegressors = zeros(ntotal,2);
        eidx = 1;
        for r = 1:length(S.functional)
            if runsToUse(r)
                nhere = size(eventRegressorsPerRun{r},1);
                erange = eidx:(eidx+nhere-1);
                eventRegressors(erange,:) = eventRegressorsPerRun{r};
                otherRegressors(erange,:) = otherRegressorsPerRun{r};
                eidx = eidx + nhere;
            end
        end
            
        %% crop to the ones we use (use last P loaded);
        
        % get regressor indices
        load('names.mat'); % sentences+stimwords
        [check,indices_stimwords] = ismember(stimwords,P.names_all);
        if sum(check==0); fprintf('error: stimwords missing regressors!\n');end
        [check,indices_sentences] = ismember(sentences,P.names_all);
        if sum(check==0); fprintf('error: sentences missing regressors!\n');end
        indices = [indices_stimwords(:);indices_sentences(:)];
        eventRegressorNames = {stimwords{:},sentences{:}};
        
        % crop
        eventRegressors = eventRegressors(:,indices);
        %otherRegressors; no need
        
        %% create data
        
        fprintf('loading imaging data\n');
        
        % load 3D mask for this subject
        %V = spm_vol(fullfile(EXPT.analysis_dir,S.name,sprintf('model3'),'mask_allruns.img'));
        V = spm_vol(fullfile(EXPT.analysis_dir,S.name,sprintf('%s/run1',M),'mask.img'));
        volmask = spm_read_vols(V);
        volmask_original = volmask;
        indicesIn3D_original = find(volmask(:));
        
        %% load the AAL atlas 
        fpath = which('atlas_aal_resized.nii');
        V = spm_vol(fpath);
        volaal = spm_read_vols(V);
        fpath = which('atlas_aal_resized.txt');
        [discard,labels_aal] = textread(fpath,'%d\t%s','delimiter','\n');
        
        %% load the language localizer group parcels
        fpath = which('langloc_parcels.nii');
        V = spm_vol(fpath);
        vollangloc = spm_read_vols(V);
        fpath = which('langloc_parcels.txt');
        [discard,labels_langloc] = textread(fpath,'%d\t%s','delimiter','\n');
        
        %% load the language localizer significance masks for this subject
        S = EXPT.subject(subj);
        fpath = fullfile(EXPT.analysis_dir,S.name,'localizers','S_N.mat');
        if exist(fpath,'file')
            % OK
        else
            fpath = fullfile(EXPT.analysis_dir,S.name,'localizers','S-N.mat');
            if exist(fpath,'file')
                % OK
            else
                fprintf('error: cannot find S-N.mat or S_N.mat in localizers\n');return;
            end
        end
        load(fpath); % gets us p, same number of values as mask

        significant_mask_0_001 = p <= 0.001;
        significant_mask_0_01  = p <= 0.01;
        
        volsig_0_001 = zeros(size(volmask));
        volsig_0_001(indicesIn3D_original) = significant_mask_0_001;
        volsig_0_01 = zeros(size(volmask));
        volsig_0_01(indicesIn3D_original) = significant_mask_0_01;
        
        %% now recalculate the subject mask (AND with AAL)
        %% and use it to remove voxels from the others
        
        volmask = volmask & (volaal > 0);
        volaal  = volaal  .* volmask;
        
        vollangloc   = vollangloc   .* volmask;
        volsig_0_001 = volsig_0_001 .* volmask;
        volsig_0_01  = volsig_0_01  .* volmask;
        
        %% and create multimasks from them
        
        indicesIn3D = find(volmask(:));
        nvoxels = length(indicesIn3D);
        
        multimask_aal   = volaal(indicesIn3D);
        multimask_group = vollangloc(indicesIn3D);
        multimask_0_001 = volsig_0_001(indicesIn3D);
        multimask_0_01  = volsig_0_01(indicesIn3D);
        
        varsToSave1 = 'multimask_aal multimask_group multimask_0_01 multimask_0_001';
        varsToSave2 = 'volmask volaal vollangloc volsig_0_01 volsig_0_001 labels_aal labels_langloc volmask_original';

        % load data
        images = zeros(ntotal,nvoxels);

        eidx = 1;
        sname = S.name;
        for r = 1:length(S.functional)
            niftidir = S.functional(r).niftidir;
            rname    = S.functional(r).name;
            run      = S.functional(r).run;
            date     = S.functional(r).date;
            session  = S.functional(r).session;
            rlength  = S.functional(r).length; % in TRs
            
            if runsToUse(r)
 
                files = fmri_get(fullfile(niftidir,sprintf('wf%s_%s-%3.4d-*.nii',sname,date,run)));
                fprintf('\trun %d\t%d files\n',r,length(files));

                for e = 1:length(files)
                    V = spm_vol(files(e,:));
                    vol = spm_read_vols(V);
                    images(eidx,:) = vol(indicesIn3D);
                    eidx = eidx + 1;
                end
                
                %nhere = size(eventRegressorsPerRun{r},1);
                %erange = eidx:(eidx+nhere-1);
                %eventRegressors(erange,:) = eventRegressorsPerRun{r};
                %otherRegressors(erange,:) = otherRegressorsPerRun{r};
                %eidx = eidx + nhere;
            end
        end
        
        % convert to percent-signal-change
        % (relative to the mean image of each session)

        fprintf('converting to PSC\n');
        for s = sessions
            fprintf('\tsession %d\n',s);
            indices = find(otherRegressors(:,2)==s);
            ninsession = length(indices);
            
            sessionMeanImage{s} = mean(images(indices,:),1);
            maskok = (sessionMeanImage{s} ~= 0);
            vin = find(maskok); vout = find(~maskok);

            tmp = repmat(sessionMeanImage{s}(vin),ninsession,1);
            
            % these voxels are OK
            images(indices,vin) = 100*((images(indices,vin) - tmp) ./ tmp);
            
            % these voxels are 0 on average, treat separately
            images(indices,vout) = 0;
            
            clear tmp;
        end
        
        %% output the PSC

        outputFile = fullfile(EXPT.analysis_dir,S.name,'imagesMPSC.mat');
        if 1
            varsToSave3 = 'eventRegressors otherRegressors eventRegressorNames images sessionMeanImage';
            eval(sprintf('save -v7.3 %s %s %s %s;',outputFile,varsToSave1,varsToSave2,varsToSave3));
        else
            load(outputFile);
        end

        %% generate examples
        
        %eventMask    = sum(eventRegressors,2);
        %tmp = diff([0;eventMask]>0);

        %indicesStart = find(tmp==1);
        
        
      case 'createDatasetCMUstyle'
        
        model = 3;
        
        fprintf('load MPSC\n');
        inputFile = fullfile(EXPT.analysis_dir,S.name,'imagesMPSC.mat');
        if 1
            load(inputFile);
            nvoxels = size(images,2);
        else
            % debug
            %nvoxels = 188298;
            nvoxels = 100; volmask = [];
            load('tmp.mat');
        end
                
        [nimages,nregressors] = size(eventRegressors);
        
        % collapse all event regressors into one
        % and determine where each trial starts
        allEvents = sum(eventRegressors,2);
        maskStart = diff([0;allEvents>0])>0; % 1 at start points
        ntrials = sum(maskStart);

        
        %% create examples
        
        examples = zeros(ntrials,nvoxels);
        labels   = zeros(ntrials,1);

        % to account for haemodynamic response
        TRshifts = [2 3];
        % how many TRs to average
        TRwindows = 2;
        
        for TRshift = TRshifts
            for TRwindow = TRwindows
                
                fprintf('creating examples TRshift=%d TRwindow=%d\n',TRshift,TRwindow);
                eidx = 1;
                for e = 1:nimages
                    if maskStart(e)
                        % trial starts here, find which regressor
                        labels(eidx) = find(eventRegressors(e,:)>0);

                        % compute range to actually use
                        if allEvents(e) <= 0.2
                            first = e+1; % assume it starts on the next
                        else
                            first = e; % starts here
                        end
                        first = first + TRshift; % shift for HDR
                        erange = first:(first+TRwindow-1);
                        
                        fprintf('\t%d\t(%1.2f)\t%d\t%d-%d\t%s\n',labels(eidx),allEvents(e),e,erange(1),erange(end),eventRegressorNames{labels(eidx)});
                        
                        % average images
                        examples(eidx,:) = mean(images(erange,:),1);
                        
                        eidx = eidx + 1;
                    end
                end

                
                %% create separate examples for stimwords and sentences

                % gives "sentences" and "stim_words"
                load('stimulus_names.mat');
                [present_sentences,indices_sentences] = ismember(eventRegressorNames,sentences);
                [present_stimwords,indices_stimwords] = ismember(eventRegressorNames,stimwords);
                nsentences = length(sentences);
                nstimwords = length(stimwords);
       
                % create examples for sentences
                examples_sentences = zeros(nsentences,nvoxels);
                for is = 1:nsentences
                    regressor = find(indices_sentences == is);
                    indices   = find(labels == regressor);
                    fprintf('\tsentence %d (%s)\nregressor %d (%s) -> %d\n',is,sentences{is},regressor,eventRegressorNames{regressor},length(indices));
            
                    examples_sentences(is,:) = mean(examples(indices,:),1);
                end
                
                % create examples for stimwords
                examples_stimwords = zeros(nstimwords,nvoxels);
                for is = 1:nstimwords
                    regressor = find(indices_stimwords == is);
                    indices   = find(labels == regressor);
                    fprintf('\tword %d (%s)\nregressor %d (%s) -> %d\n',is,stimwords{is},regressor,eventRegressorNames{regressor},length(indices));
            
                    examples_stimwords(is,:) = mean(examples(indices,:),1);
                end 
        
        
                %% output

                % these come from the MPSC file
                varsToSave1 = 'multimask_aal multimask_group multimask_0_01 multimask_0_001';
                varsToSave2 = 'volmask volaal vollangloc volsig_0_01 volsig_0_001 labels_aal labels_langloc volmask_original';
                
                outputFile = fullfile(EXPT.analysis_dir,S.name,sprintf('examplesMPSC_S%d_A%d.mat',TRshift,TRwindow));
                if 1
                    varsToSave3 = 'eventRegressorNames volmask TRshift TRwindow examples labels examples_sentences sentences examples_stimwords stimwords';
                    eval(sprintf('save -v7.3 %s %s %s %s;',outputFile,varsToSave1,varsToSave2,varsToSave3));
                else
                    load(outputFile);
                end
                
            end; % loop over TR windows
        end; % loop over TR shift
            
    end; % end of switch