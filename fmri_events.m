function names = fmri_events(EXPT,model,subj,runname,getall)
    
    % Get event names.
    %
    % USAGE: names = fmri_events(EXPT,model,subj,runname,[getall])
    %
    % INPUTS:
    %   EXPT - experiment structure
    %   model - model number
    %   subj - subject number
    %   runname - name of run
    %   getall (optional) - if true, will get all event names in the para
    %           file; if false, will only get names of events that occur in the run
    %           (deafult: true)
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
    
    if nargin < 4 || getall
        names = P.names_all';
    else
        names = P.names';
    end