function names = parse_aal
    
    fdir = fileparts(which('aal_MNI_V4.txt'));
    T = fullfile(fdir,'aal_MNI_V4.txt');
    f = fopen(T);
    fgetl(f);   % remove header
    while ~feof(f)
        L = fgetl(f);
        words = regexp(L,'\t','split');
        n = str2num(words{1});
        names{n,1} = words{2};
    end
    
    fclose(f);