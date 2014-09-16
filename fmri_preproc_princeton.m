function fmri_preproc(EXPT,subj,tasks)
    
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
                disp(['run ',num2str(r)]);
                dicomdir = S.functional(r).dicomdir;
                niftidir = S.functional(r).niftidir;
                run = S.functional(r).run;
                files = dir(fullfile(dicomdir,sprintf('*-%d-*',run)));
                files = dir2char(files,dicomdir);
                hdr = spm_dicom_headers(files);
                if ~exist(niftidir,'dir'); mkdir(niftidir);  end
                delete(fullfile(niftidir,'w*'));
                delete(fullfile(niftidir,'s*'));
                cd(niftidir);
                spm_dicom_convert(hdr,'all','flat','nii');
            end
            disp('anatomical');
            dicomdir = S.anatomical.dicomdir;
            niftidir = S.anatomical.niftidir;
            run = S.anatomical.run;
            files = dir(fullfile(dicomdir,sprintf('*-%d-*',run)));
            files = dir2char(files,dicomdir);
            hdr = spm_dicom_headers(files);
            if ~exist(niftidir,'dir'); mkdir(niftidir); end
            cd(niftidir);
            spm_dicom_convert(hdr,'all','flat','nii');
            cd(curdir);
        
        case 'realign'
            % This step does motion correction. Only a mean image (mean_*) is
            % resliced; all other images have their headers modified.
            
            disp('Realigning...')
            
            % get nifti filenames
            for r = 1:length(S.functional)
                niftidir = S.functional(r).niftidir;
                run = S.functional(r).run;
                files{r,1} = fmri_get(fullfile(niftidir,sprintf('f*-%3.4d-*',run)));
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

        case 'coregister'
            % This step first coregisters the mean functional to the
            % anatomical, and then coregisters the anatomical to the
            % MNI template.
            
            disp('Coregistration...')
            
            % 1) mean functional -> anatomical
            % use normalied mutual information for registering images from
            % different modalities (epi -> T1)
            niftidir = S.functional(1).niftidir;
            mean_epi = fmri_get(fullfile(niftidir,'mean*'));
            run = S.anatomical.run;
            anatomical = fmri_get(fullfile(S.anatomical.niftidir,sprintf('*-%3.4d-*',run)));
            T1 = spm_coreg(anatomical,mean_epi,struct('cost_fun','nmi','graphics',0));
            
            % transform all other functionals
            for r = 1:length(S.functional)
                niftidir = S.functional(r).niftidir;
                run = S.functional(r).run;
                P{r,1} = fmri_get(fullfile(niftidir,sprintf('*-%3.4d-*',run)));
            end
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
            
        case 'normalize'
            % Segment and normalize anatomical image to MNI template.
            % Functionals and anatomical prepended with 'w'.
            
            disp('Normalizing anatomical...')
            
            run = S.anatomical.run;
            anatomical = fmri_get(fullfile(S.anatomical.niftidir,sprintf('*-%3.4d-*',run)));
            res = spm_preproc(anatomical);   % compute warping parameters
            sn = spm_prep2sn(res);
            spm_write_sn(anatomical,sn);     % normalize anatomical using the warp parameters already calculated
            
            % normalize functionals
            for r = 1:length(S.functional)
                niftidir = S.functional(r).niftidir;
                run = S.functional(r).run;
                P = fmri_get(fullfile(niftidir,sprintf('f*-%3.4d-*',run)));
                spm_write_sn(P,sn);
            end
            
            adir = fullfile(EXPT.analysis_dir,EXPT.subject(subj).name);
            if ~exist(adir,'dir'); mkdir(adir); end
            save(fullfile(adir,'normalization_params'),'sn');
            
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
    end