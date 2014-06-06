function names = fmri_events(EXPT,model,subj,runname)
    
    % Get event names.
    %
    % USAGE: names = fmri_events(EXPT,model,subj,runname)
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number
    %   subj - subject number
    %   runname - name of run
    %
    % OUTPUTS:
    %   names - cell array of event names
    %
    % Sam Gershman, June 2014
    
    S = EXPT.subject(subj);
    
    for run = 1:length(S.functional)
        if strcmp(S.functional(run).name,runname); break; end
    end
    
    P = parse_para(S.functional(run).para{model},EXPT.TR);
    names = P.names_all';