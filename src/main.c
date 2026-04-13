#include "solver.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static SolverConfig default_config(void) {
    SolverConfig cfg;
    cfg.tol_x = 1e-10;
    cfg.tol_f = 1e-10;
    cfg.max_iter = 300;
    cfg.scan_steps = 3000;
    return cfg;
}

static int parse_coefficients_line(const char *line, SexticCoefficients *coeffs) {
    double values[7] = {0};
    char *end = NULL;
    const char *cursor = line;

    for (int i = 0; i < 7; ++i) {
        while (isspace((unsigned char)*cursor)) {
            ++cursor;
        }

        if (*cursor == '\0' || *cursor == '\n') {
            return 0;
        }

        values[i] = strtod(cursor, &end);
        if (end == cursor) {
            return 0;
        }
        cursor = end;
    }

    while (isspace((unsigned char)*cursor)) {
        ++cursor;
    }

    if (*cursor != '\0' && *cursor != '\n') {
        return 0;
    }

    coeffs->a = values[0];
    coeffs->b = values[1];
    coeffs->c = values[2];
    coeffs->d = values[3];
    coeffs->n = values[4];
    coeffs->m = values[5];
    coeffs->k = values[6];
    return 1;
}

static void read_coefficients_interactive(SexticCoefficients *coeffs) {
    char line[512];

    for (;;) {
        printf("Enter coefficients for ax^6 + bx^5 + cx^4 + dx^3 + nx^2 + mx + k = 0\n");
        printf("Format: a b c d n m k\n> ");

        if (fgets(line, sizeof(line), stdin) == NULL) {
            fprintf(stderr, "Input stream ended.\n");
            exit(EXIT_FAILURE);
        }

        if (parse_coefficients_line(line, coeffs)) {
            return;
        }

        printf("Invalid input. Please enter exactly 7 numeric values.\n\n");
    }
}

static void print_roots(const RealRoots *roots) {
    printf("Roots: [");
    for (size_t i = 0; i < roots->count; ++i) {
        if (i > 0) {
            printf(", ");
        }
        double iv = (double)((long long)roots->values[i]);
        if (roots->values[i] == iv) {
            printf("%.0f", roots->values[i]);
        } else {
            printf("%.6f", roots->values[i]);
        }
    }
    printf("]\n");
}

static int ask_continue(void) {
    char line[32];
    for (;;) {
        printf("Continue? (y/n): ");
        if (fgets(line, sizeof(line), stdin) == NULL) {
            return 0;
        }
        if (line[0] == 'y' || line[0] == 'Y') {
            return 1;
        }
        if (line[0] == 'n' || line[0] == 'N') {
            return 0;
        }
        printf("Please type 'y' or 'n'.\n");
    }
}

int main(void) {
    SolverConfig cfg = default_config();

    for (;;) {
        SexticCoefficients coeffs;
        RealRoots roots;
        char err_buf[256] = {0};

        read_coefficients_interactive(&coeffs);

        SolveStatus status = solve_equation(&coeffs, &cfg, &roots, err_buf, sizeof(err_buf));
        if (status == SOLVE_ERROR) {
            printf("Error: %s\n", err_buf);
        } else if (status == SOLVE_ANY_NUMBER) {
            if (coeffs.a == 0.0 && coeffs.b == 0.0 && coeffs.c == 0.0 && coeffs.d == 0.0 && coeffs.n == 0.0 &&
                coeffs.m == 0.0 && coeffs.k == 0.0) {
                printf("Any number is a solution.\n");
            } else {
                printf("Any number.\n");
            }
        } else if (status == SOLVE_NO_SOLUTIONS) {
            printf("No solutions.\n");
        } else {
            print_roots(&roots);
        }

        if (!ask_continue()) {
            break;
        }
        printf("\n");
    }

    return EXIT_SUCCESS;
}
