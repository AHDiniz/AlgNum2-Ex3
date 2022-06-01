function iterations = calc_iter(iter, n)
    [nr, nc] = size(iter);
    if nr == 1 && nc == 1
        iterations = iter;
        return;
    endif
    if nr == 0 && nc == 0
        iterations = 0;
        return;
    endif
    iterations = (iter(1, 1) - 1) * n + iter(1, 2);
endfunction