function [c n] = find_regressors(regnames, name)
    
    % Get a regressor mask of those matching name.
    %
    % USAGE: [c n] = find_regressors(regnames, name)
    %
    % INPUTS:
    %   regnames - cell array of regressor names
    %   name - string or cell array of names
    %
    % OUTPUTS:
    %   c - regressor mask
    %   n - mask identifier
    
    if isempty(name)
        c = zeros(1, length(regnames));
        n = '';
    elseif isnumeric(name) && isvector(name) && length(name) == length(regnames)
        c = name;
        n = [horzcat(regnames{c>0}) '>' horzcat(regnames{c<0})];
    elseif iscell(name)
        [c nt] = cellfun(@(nt) find_regressors(regnames, nt), name(:), 'UniformOutput', 0);
        c = any(cell2mat(c));
        n = nt{1};
        for i = 2:length(nt)
            n = [n '+' nt{i}];
        end
    elseif isstruct(name)
        n = name.name;
        c = name.c;
    else
        n = name;
        if strncmp(n, '^Sn\(\d+\) ', 11)
            n = n(12:end);
        end
        if length(n) >= 1 && n(end) == '$'
            n = n(1:end-1);
        end
        if length(n) >= 9 && strcmp(n(end-8:end), '\*bf\(1\)')
            n = n(1:end-9);
        end
        if length(n) >= 3 && strcmp(n(end-2:end), '\^1')
            n = n(1:end-3);
        end
        if name(1) ~= '^'
            name = ['^Sn\(\d+\) (.*x)?(' name ')(\^1)?(\*bf\(1\))?$'];
        end
        c = ~cellfun('isempty',regexp(regnames, name));
        if ~any(c)
            warning('find_regressors:notfound', 'No matching regressors for "%s"', name);
        end
    end
