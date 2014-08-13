function allRegressors = get_all_regressors(subj,model)
    % Get cell array of all regressor names for a subject/model.
    %
    % USAGE: allRegressors = get_all_regressors(subj,model)
    %
    % INPUTS:
    %   subj - subject number
    %   model - model number
    %
    % OUTPUTS:
    %   allRegressors - cell array of regressor names
    %
    % Walid Bendris, June 2014
	
    curDir=cd;
    cd(['/srv/old_data/projects/MACKEREL/analysis02/subj' sprintf('%02d', subj) '/model' num2str(model)]);
    load('SPM.mat');
    allRegressors={}; jStart=1;
    len=str2num(SPM.xX.name{end}(4:end-10));
    for i=1:length(SPM.xX.name)
        for j=jStart:len
            if strncmp(SPM.xX.name{i}, ['Sn(' num2str(j) ') '], 6) && ...
                    strcmp(SPM.xX.name{i}(end-5:end), '*bf(1)')
                    allRegressors{end+1,1}= SPM.xX.name{i}(7:end-6);
                    if j>jStart, jStart=j; end
                    break;
            elseif strncmp(SPM.xX.name{i}, ['Sn(' num2str(j) ') '], 7) && ...
                    strcmp(SPM.xX.name{i}(end-5:end), '*bf(1)')
                    allRegressors{end+1,1}= SPM.xX.name{i}(8:end-6);
                    if j>jStart, jStart=j; end
                    break;
            end
        end
    end, allRegressors=unique(allRegressors);
    cd curDir;
end