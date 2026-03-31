#include "solver.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static bool is_finite_number(double x) {
    return isfinite(x);
}

bool validate_coefficients(const SexticCoefficients *coeffs, char *err_buf, size_t err_size) {
    if (coeffs == NULL) {
        snprintf(err_buf, err_size, "Internal error: null coefficient pointer.");
        return false;
    }

    const double values[] = {coeffs->a, coeffs->b, coeffs->c, coeffs->d, coeffs->n, coeffs->m, coeffs->k};
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); ++i) {
        if (!is_finite_number(values[i])) {
            snprintf(err_buf, err_size, "Coefficient %zu is not a finite real number.", i + 1);
            return false;
        }
    }

    if (fabs(coeffs->a) < 1e-15) {
        snprintf(err_buf, err_size,
                 "Coefficient 'a' is zero (or too close). This is not a sixth-degree equation.");
        return false;
    }

    return true;
}

double eval_polynomial(const SexticCoefficients *coeffs, double x) {
    return ((((((coeffs->a * x + coeffs->b) * x + coeffs->c) * x + coeffs->d) * x + coeffs->n) * x + coeffs->m) * x +
            coeffs->k);
}

double eval_derivative(const SexticCoefficients *coeffs, double x) {
    return (((((6.0 * coeffs->a * x + 5.0 * coeffs->b) * x + 4.0 * coeffs->c) * x + 3.0 * coeffs->d) * x +
             2.0 * coeffs->n) *
                x +
            coeffs->m);
}

double cauchy_root_bound(const SexticCoefficients *coeffs) {
    const double a_abs = fabs(coeffs->a);
    const double others[] = {fabs(coeffs->b), fabs(coeffs->c), fabs(coeffs->d), fabs(coeffs->n), fabs(coeffs->m),
                             fabs(coeffs->k)};

    double max_other = 0.0;
    for (size_t i = 0; i < sizeof(others) / sizeof(others[0]); ++i) {
        if (others[i] > max_other) {
            max_other = others[i];
        }
    }

    return 1.0 + max_other / a_abs;
}

bool find_real_roots(const SexticCoefficients *coeffs, const SolverConfig *cfg, RealRoots *out_roots,
                     char *err_buf, size_t err_size) {
    (void)coeffs;
    (void)cfg;

    if (out_roots == NULL) {
        snprintf(err_buf, err_size, "Internal error: null roots output pointer.");
        return false;
    }

    out_roots->count = 0;

    /*
     * Project skeleton note:
     * TODO: Implement robust interval isolation + bisection/Newton hybrid.
     * TODO: Add duplicate-root suppression with tolerance.
     */
    snprintf(err_buf, err_size, "Root-finding engine not implemented yet (project skeleton).");
    return false;
}
