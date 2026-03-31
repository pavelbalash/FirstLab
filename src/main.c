#include "solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static SolverConfig default_config(void) {
    SolverConfig cfg;
    cfg.tol_x = 1e-10;
    cfg.tol_f = 1e-10;
    cfg.max_iter = 200;
    cfg.scan_steps = 2000;
    return cfg;
}

static int read_coefficients(SexticCoefficients *coeffs) {
    printf("Enter coefficients for ax^6 + bx^5 + cx^4 + dx^3 + nx^2 + mx + k = 0\n");
    printf("Format: a b c d n m k\n> ");

    if (scanf("%lf %lf %lf %lf %lf %lf %lf", &coeffs->a, &coeffs->b, &coeffs->c, &coeffs->d, &coeffs->n,
              &coeffs->m, &coeffs->k) != 7) {
        return 0;
    }
    return 1;
}

static void print_demo_example(void) {
    SexticCoefficients coeffs = {.a = 1.0, .b = -21.0, .c = 175.0, .d = -735.0, .n = 1624.0, .m = -1764.0, .k = 720.0};

    printf("Demo polynomial:\n");
    printf("x^6 - 21x^5 + 175x^4 - 735x^3 + 1624x^2 - 1764x + 720 = 0\n");
    printf("(Expected roots: 1,2,3,4,5,6)\n\n");

    printf("Cauchy bound estimate: |x| <= %.6f\n", cauchy_root_bound(&coeffs));
}

int main(int argc, char **argv) {
    if (argc > 1 && strcmp(argv[1], "--demo") == 0) {
        print_demo_example();
        return EXIT_SUCCESS;
    }

    SexticCoefficients coeffs;
    SolverConfig cfg = default_config();
    RealRoots roots;
    char err_buf[256] = {0};

    if (!read_coefficients(&coeffs)) {
        fprintf(stderr, "Input error: expected 7 real numbers.\n");
        return EXIT_FAILURE;
    }

    if (!validate_coefficients(&coeffs, err_buf, sizeof(err_buf))) {
        fprintf(stderr, "Validation error: %s\n", err_buf);
        return EXIT_FAILURE;
    }

    printf("Search interval estimate (Cauchy): [%.6f, %.6f]\n", -cauchy_root_bound(&coeffs), cauchy_root_bound(&coeffs));

    if (!find_real_roots(&coeffs, &cfg, &roots, err_buf, sizeof(err_buf))) {
        fprintf(stderr, "Solver status: %s\n", err_buf);
        fprintf(stderr, "Tip: implement find_real_roots() in src/solver.c\n");
        return EXIT_FAILURE;
    }

    printf("Found %zu real root(s):\n", roots.count);
    for (size_t i = 0; i < roots.count; ++i) {
        printf("  root[%zu] = %.12f\n", i, roots.values[i]);
    }

    return EXIT_SUCCESS;
}
