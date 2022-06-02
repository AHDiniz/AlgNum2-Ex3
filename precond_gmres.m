infiles = {"in/cavity05.mat", "in/cz2548.mat", "in/epb3.mat"};
names = {"cavity05", "cz2548", "epb3"};
ks = {200, 200, -1};

for i = 1 : 3
    load(infiles{i});
    name = names{i};
    A = Problem.A;
    n = rows(A);
    b = ones(n, 1);

    zeroL = []; zeroU = [];
    zeroLr = []; zeroUr = [];
    croutL = []; croutU = [];
    croutLr = []; croutUr = [];

    zero_worked = zero_rcm_worked = crout_worked = crout_rcm_worked = 1;

    % No preconditioning:
    nopcd = [[0], -1, inf, [0], [0]];
    nopcd_t = 0;
    if ks{i} == -1
        k_options = [2, 3, 5, 10, 50, 100, 150, 200];
        nopcd = test_gmres(A, b, k_options, 10e-11, 1000);
        ks{i} = nopcd(6);
        nopcd_t = nopcd(7);
    else
        tic();
        nopcd = gmres(A, b, ks{i}, 10e-11, 1000);
        nopcd_t = toc();
    endif

    % ILU Zero:
    zero = [[0], -1, inf, [0], [0]];
    zero_t = 0;
    try
        ilu_opts.type = "nofill";
        [zL, zU] = ilu(A, ilu_opts);
        tic();
        zero = gmres(A, b, ks{i}, 10e-11, 1000, zL, zU);
        zero_t = toc();
        zeroL = zL;
        zeroU = zU;
    catch
        zero_worked = 0;
    end_try_catch

    % ILU Zero RCM:
    zero_rcm = [[0], -1, inf, [0], [0]];
    zero_rcm_t = 0;
    try
        perm = symrcm(A);
        I = speye(A);
        P = I(perm,:);
        R = P * A * P';
        ilu_opts.type = "nofill";
        [zLr, zUr] = ilu(R, ilu_opts);
        tic();
        zero_rcm = gmres(R, b, ks{i}, 10e-11, 1000, zLr, zUr);
        zero_rcm_t = toc();
        zeroLr = zLr;
        zeroUr = zUr;
    catch
        zero_rcm_worked = 0;
    end_try_catch

    % Crout:
    crout = [[0], -1, inf, [0], [0]];
    crout_t = 0;
    try
        crout_opts.type = "crout";
        crout_opts.droptol = 10e-4;
        [crL, crU] = ilu(A, crout_opts);
        tic();
        crout = gmres(A, b, ks{i}, 10e-11, 1000, crL, crU);
        crout_t = toc();
        croutL = crL;
        croutU = crU;
    catch
        crout_worked = 0;
    end_try_catch

    % Crout RCM:
    crout_rcm = [[0], -1, inf, [0], [0]];
    crout_rcm_t = 0;
    try
        crout_opts.type = "crout";
        crout_opts.droptol = 10e-4;
        perm = symrcm(A);
        I = speye(A);
        P = I(perm,:);
        R = P * A * P';
        [crLr, crUr] = ilu(R, crout_opts);
        tic();
        crout_rcm = gmres(A, b, ks{i}, 10e-11, 1000, crLr, crUr);
        crout_rcm_t = toc();
        croutLr = crLr;
        croutUr = crUr;
    catch
        crout_rcm_worked = 0;
    end_try_catch

    % Text outputs:
    outputs = {
        {"NO PCD", nopcd, nopcd_t},
        {"ILU Zero", zero, zero_t},
        {"ILU Zero RCM", zero_rcm, zero_rcm_t},
        {"ILU Crout", crout, crout_t},
        {"ILU Crout RCM", crout_rcm, crout_rcm_t}
    };
    print_text_outputs(name, A, n, outputs);

    % Graphics:
    fig = figure();
    plot(calc_iter(nopcd(4)), log(nopcd(5)));
    if zero_worked
        plot(calc_iter(zero(4)), log(zero(5)));
    endif
    if zero_rcm_worked
        plot(calc_iter(zero_rcm(4)), log(zero_rcm(5)));
    endif
    if crout_worked
        plot(calc_iter(crout(4)), log(crout(5)));
    endif
    if crout_rcm_worked
        plot(calc_iter(crout_rcm(4)), log(crout_rcm(5)));
    endif
    title(name);
    xlabel("Iterações");
    ylabel("log(Resíduo Relativo)");
    print(fig, strcat("out/", name, "_relres_iter.svg"), "-dsvg");

endfor
