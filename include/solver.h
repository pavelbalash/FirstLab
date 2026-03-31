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

bool validate_coefficients(const SexticCoefficients *coeffs, char *err_buf, size_t err_size);
double eval_polynomial(const SexticCoefficients *coeffs, double x);
double eval_derivative(const SexticCoefficients *coeffs, double x);
double cauchy_root_bound(const SexticCoefficients *coeffs);

bool find_real_roots(const SexticCoefficients *coeffs, const SolverConfig *cfg, RealRoots *out_roots,
                     char *err_buf, size_t err_size);

#endif
