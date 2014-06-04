function P = parse_para(filename,TR)
    
    % Parse *.para files, which contain information about onsets, durations
    % and names of events in a functional run.
    %
    % USAGE: P = parse_para(filename,[TR])
    %
    % INPUTS:
    %   filename - name of *.para file
    %   TR (optional) - repetition time; if not supplied, it is assumed
    %   that time is measured in seconds
    %
    % OUTPUTS:
    %   P - structure with the following fields:
    %           .onsets - cell array, where each cell contains the vector of
    %           onsets (in seconds) for the corresponding event
    %           .names - cell array containing the name of each event
    %           .durations - vector containing the duration (in seconds) of each event
    %           .events - vector containing event indices
    %
    % NOTE:
    %   *.para files have the following format: there are three segments of
    %   the file (onsets, durations, names), flagged by "#" (e.g.,
    %   #onsets). In the onsets segment, each line contains two numbers.
    %   The first number is the onset (either in seconds or TRs) and the
    %   second number is the event index. The names segment consists of a
    %   single line, containing the name of each event separated by spaces.
    %   The durations segment also consists of a single line, containing
    %   the duration (in seconds or TRs) of each event.
    % 
    % Sam Gershman, Jan 2014
    
    f = fopen(filename);
    ons = [];
    event = [];
    
    while ~feof(f)
        text = fgetl(f);
        text = strtrim(text);
        if ~isempty(text)
            if text(1) == '#'
                fieldname = text(2:end);
            else
                w = regexp(text,' ','split');
                switch fieldname
                    case 'onsets'
                        ons(end+1) = str2num(w{1});
                        event(end+1) = str2num(w{2});
                    case 'names'
                        P.names = w;
                    case 'durations'
                        for i = 1:length(w)
                            P.durations(i) = str2num(w{i});
                        end
                end
            end
        end
    end
    
    fclose(f);
    
    % convert to seconds if necessary
    if nargin > 1
        ons = ons*TR;
        P.durations = P.durations*TR;
    end
    
    u = unique(event);
    P.events = u;
    P.onsets = cell(1,length(u));
    for i = 1:length(u)
        P.onsets{i} = ons(event==u(i));
    end
    P.names_all = P.names;
    P.names = P.names(u);
    P.durations = P.durations(u);