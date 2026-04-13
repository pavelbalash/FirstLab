#ifndef SOLVER_H
#define SOLVER_H

#include <stdbool.h>
#include <stddef.h>

typedef struct {
    double a;
    double b;
    double c;
    double d;
    double n;
    double m;
    double k;
} SexticCoefficients;

typedef struct {
    double values[6];
    size_t count;
} RealRoots;

typedef struct {
    double tol_x;
    double tol_f;
    int max_iter;
    int scan_steps;
} SolverConfig;

typedef enum {
    SOLVE_OK = 0,
    SOLVE_ANY_NUMBER,
    SOLVE_NO_SOLUTIONS,
    SOLVE_ERROR
} SolveStatus;

bool validate_coefficients(const SexticCoefficients *coeffs, char *err_buf, size_t err_size);
SolveStatus solve_equation(const SexticCoefficients *coeffs,
                           const SolverConfig *cfg,
                           RealRoots *out_roots,
                           char *err_buf,
                           size_t err_size);

#endif
