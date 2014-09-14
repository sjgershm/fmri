function names = parse_aal
    
    spmdir = fileparts(which('SPM.m'));
    T = fullfile(spmdir,'toolbox/WFU_PickAtlas_3.0.4/wfu_pickatlas/MNI_atlas_templates/aal_MNI_V4.txt');
    f = fopen(T);
    fgetl(f);   % remove header
    while ~feof(f)
        L = fgetl(f);
        words = regexp(L,'\t','split');
        n = str2num(words{1});
        names{n,1} = words{2};
    end
    
    fclose(f);