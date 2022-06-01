function [x, f, relres, iter, resvec, k, elapsed] = test_gmres(A, b, k_options, tol, max_iter)
    best_relres = inf;
    for i = 1 : numel(k_options)
        tic();
        [xi, fi, relresi, iteri, resveci] = gmres(A, b, k_options(i), tol, max_iter);
        t = toc();
        if relresi < best_relres
            x = xi;
            f = fi;
            relres = relresi;
            iter = iteri;
            resvec = resveci;
            k = k_options(i);
            elapsed = t;
        endif
    endfor
endfunction