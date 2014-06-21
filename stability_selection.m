function [rank, c] = stability_selection(Y,rep)
    
    V = size(Y,2);
    c = zeros(1,V);
    ii = 0;
    u = unique(rep);
    N = length(u);
    for n = 1:N
        y = Y(rep==u(n),:);
        dn = bsxfun(@minus,y,mean(y));
        fn = sum(dn.^2);
        for m = (n+1):N
            ii = ii + 1;
            y = Y(rep==u(m),:);
            dm = bsxfun(@minus,y,mean(y));
            fm = sum(dm.^2);
            c = c + sum(dn.*dm)./sqrt(fm.*fn);
        end
    end
    
    [c,rank] = sort(c/ii,'descend');