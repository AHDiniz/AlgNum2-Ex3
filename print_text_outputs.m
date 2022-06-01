function print_text_outputs(filename, A, n, outputs)
    f = fopen(strcat("out/", filename, ".txt"), "w");
    fprintf(f, "Matriz: %s\n", filename);
    fprintf(f, "Ordem do sistema = %d\n", n);
    fprintf(f, "Não nulos = %d\n\n", nnz(A));
    for i = 1 : numel(outputs)
        name = outputs{i}{1};
        printf("%s\n", name);
        result = outputs{i}{2};
        printf("%d\n", columns(result));
        t = outputs{i}{3};
        fprintf(f, "%s:\n", name);
        fprintf(f, "- Norma da solução: %f\n", norm(result(1), 2));
        fprintf(f, "- Resíduo relativo: %f\n", result(3));
        fprintf(f, "- Número de iterações: %d\n", calc_iter(result(4), n));
        fprintf(f, "- Tempo de execução: %fs\n", t);
        fprintf(f, "- Flag de convergência: %d\n\n", result(2));
    endfor
    fclose(f);
endfunction