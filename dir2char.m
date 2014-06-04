function x = dir2char(f,d)
    
    % Convert directory structure (from dir.m) into a character array.
    %
    % USAGE: x = dir2char(f,d)
    %
    % INPUTS:
    %   f - directory structure
    %   d - file path
    %
    % OUTPUTS:
    %   x - character array of file names
    %
    % Sam Gershman, Jan 2014
    
    for i = 1:length(f)
        x{i} = fullfile(d,f(i).name);
    end
    x = char(x);