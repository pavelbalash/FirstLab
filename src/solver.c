#include "solver.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#define MAX_DEGREE 6
#define EPS_ZERO 1e-12
#define PI 3.14159265358979323846

typedef struct {
    double coeffs[MAX_DEGREE + 1];
    int degree;
} Polynomial;

static double abs_max(const double *arr, size_t n) {
    double m = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double v = fabs(arr[i]);
        if (v > m) {
            m = v;
        }
    }
    return m;
}

bool validate_coefficients(const SexticCoefficients *coeffs, char *err_buf, size_t err_size) {
    if (coeffs == NULL) {
        snprintf(err_buf, err_size, "Internal error: null coefficient pointer.");
        return false;
    }

    const double values[] = {coeffs->a, coeffs->b, coeffs->c, coeffs->d, coeffs->n, coeffs->m, coeffs->k};
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); ++i) {
        if (!isfinite(values[i])) {
            snprintf(err_buf, err_size, "Coefficient %zu is not a finite real number.", i + 1);
            return false;
        }
    }
    return true;
}

static double eval_poly(const Polynomial *p, double x) {
    double acc = p->coeffs[0];
    for (int i = 1; i <= p->degree; ++i) {
        acc = acc * x + p->coeffs[i];
    }
    return acc;
}

static double eval_poly_derivative(const Polynomial *p, double x) {
    if (p->degree <= 0) {
        return 0.0;
    }

    double acc = p->coeffs[0] * p->degree;
    for (int i = 1; i < p->degree; ++i) {
        acc = acc * x + p->coeffs[i] * (p->degree - i);
    }
    return acc;
}

static void remove_leading_zeros(Polynomial *p, double eps) {
    while (p->degree > 0 && fabs(p->coeffs[0]) <= eps) {
        for (int i = 0; i < p->degree; ++i) {
            p->coeffs[i] = p->coeffs[i + 1];
        }
        p->degree -= 1;
    }
}

static double root_bound(const Polynomial *p) {
    if (p->degree <= 0) {
        return 1.0;
    }
    double lead = fabs(p->coeffs[0]);
    if (lead < EPS_ZERO) {
        return 10.0;
    }
    double m = abs_max(&p->coeffs[1], (size_t)p->degree);
    return 1.0 + m / lead;
}

static bool add_root(RealRoots *roots, double value) {
    if (roots->count >= 6) {
        return false;
    }
    roots->values[roots->count++] = value;
    return true;
}

static bool linear_method(const Polynomial *p, double *root) {
    if (p->degree != 1) {
        return false;
    }
    double a = p->coeffs[0];
    if (fabs(a) < EPS_ZERO) {
        return false;
    }
    *root = -p->coeffs[1] / a;
    return isfinite(*root);
}

static size_t vieta_method(const Polynomial *p, double *roots) {
    if (p->degree != 2) {
        return 0;
    }

    double a = p->coeffs[0];
    double b = p->coeffs[1];
    double c = p->coeffs[2];
    if (fabs(a) < EPS_ZERO) {
        return 0;
    }

    double d = b * b - 4.0 * a * c;
    if (d < -EPS_ZERO) {
        return 0;
    }
    if (fabs(d) <= EPS_ZERO) {
        roots[0] = -b / (2.0 * a);
        return 1;
    }

    double sd = sqrt(d);
    roots[0] = (-b - sd) / (2.0 * a);
    roots[1] = (-b + sd) / (2.0 * a);
    return 2;
}

static size_t cardano_method(const Polynomial *p, double *roots) {
    if (p->degree != 3) {
        return 0;
    }

    double a = p->coeffs[0], b = p->coeffs[1], c = p->coeffs[2], d = p->coeffs[3];
    if (fabs(a) < EPS_ZERO) {
        return 0;
    }

    double p_val = (3.0 * a * c - b * b) / (3.0 * a * a);
    double q_val = (2.0 * b * b * b - 9.0 * a * b * c + 27.0 * a * a * d) / (27.0 * a * a * a);
    double disc = q_val * q_val / 4.0 + p_val * p_val * p_val / 27.0;
    double shift = -b / (3.0 * a);

    if (disc > EPS_ZERO) {
        double u = cbrt(-q_val / 2.0 + sqrt(disc));
        double v = cbrt(-q_val / 2.0 - sqrt(disc));
        roots[0] = u + v + shift;
        return 1;
    }

    if (fabs(disc) <= EPS_ZERO) {
        double u = cbrt(-q_val / 2.0);
        roots[0] = 2.0 * u + shift;
        roots[1] = -u + shift;
        return 2;
    }

    double r = sqrt(-p_val * p_val * p_val / 27.0);
    if (r < EPS_ZERO) {
        roots[0] = shift;
        return 1;
    }

    double phi = acos(-q_val / (2.0 * r));
    double t = 2.0 * cbrt(r);
    roots[0] = t * cos(phi / 3.0) + shift;
    roots[1] = t * cos((phi + 2.0 * PI) / 3.0) + shift;
    roots[2] = t * cos((phi + 4.0 * PI) / 3.0) + shift;
    return 3;
}

static bool newton_method(const Polynomial *p, const SolverConfig *cfg, double *root) {
    double bound = root_bound(p);
    const int seeds = 24;

    for (int s = 0; s <= seeds; ++s) {
        double x = -bound + (2.0 * bound * s) / seeds;
        for (int i = 0; i < cfg->max_iter; ++i) {
            double fx = eval_poly(p, x);
            double dfx = eval_poly_derivative(p, x);
            if (fabs(dfx) < 1e-14) {
                break;
            }
            double nx = x - fx / dfx;
            if (!isfinite(nx)) {
                break;
            }
            if (fabs(nx - x) < cfg->tol_x && fabs(eval_poly(p, nx)) < cfg->tol_f * 50.0) {
                *root = nx;
                return true;
            }
            x = nx;
        }
    }
    return false;
}

static bool bracket_sign_change(const Polynomial *p, const SolverConfig *cfg, double *left, double *right) {
    double bound = root_bound(p);
    int steps = cfg->scan_steps > 100 ? cfg->scan_steps : 100;
    double dx = (2.0 * bound) / steps;

    double x0 = -bound;
    double f0 = eval_poly(p, x0);
    for (int i = 1; i <= steps; ++i) {
        double x1 = -bound + i * dx;
        double f1 = eval_poly(p, x1);
        if (fabs(f0) < cfg->tol_f) {
            *left = x0 - dx;
            *right = x0 + dx;
            return true;
        }
        if (f0 * f1 <= 0.0) {
            *left = x0;
            *right = x1;
            return true;
        }
        x0 = x1;
        f0 = f1;
    }
    return false;
}

static bool chord_tangent_method(const Polynomial *p, const SolverConfig *cfg, double *root) {
    double a, b;
    if (!bracket_sign_change(p, cfg, &a, &b)) {
        return false;
    }

    double fa = eval_poly(p, a);
    double fb = eval_poly(p, b);
    for (int i = 0; i < cfg->max_iter; ++i) {
        if (fabs(fa) < cfg->tol_f) {
            *root = a;
            return true;
        }
        if (fabs(fb) < cfg->tol_f) {
            *root = b;
            return true;
        }

        /* chord (secant on bracket endpoints) */
        double x_chord = (fabs(fb - fa) > 1e-14) ? (a - fa * (b - a) / (fb - fa)) : (a + b) * 0.5;
        /* tangent (Newton from the better endpoint) */
        double x_tangent = (fabs(fa) < fabs(fb)) ? a : b;
        double dft = eval_poly_derivative(p, x_tangent);
        if (fabs(dft) > 1e-14) {
            x_tangent -= eval_poly(p, x_tangent) / dft;
        }

        if (!isfinite(x_chord) || x_chord <= a || x_chord >= b) {
            x_chord = (a + b) * 0.5;
        }
        if (!isfinite(x_tangent) || x_tangent <= a || x_tangent >= b) {
            x_tangent = (a + b) * 0.5;
        }

        double x_mid = 0.5 * (x_chord + x_tangent);
        double f_mid = eval_poly(p, x_mid);

        if (fabs(f_mid) < cfg->tol_f) {
            *root = x_mid;
            return true;
        }
        if (fa * f_mid <= 0.0) {
            b = x_mid;
            fb = f_mid;
        } else {
            a = x_mid;
            fa = f_mid;
        }

        if (fabs(b - a) < cfg->tol_x) {
            *root = 0.5 * (a + b);
            return true;
        }
    }
    return false;
}

static bool secant_method(const Polynomial *p, const SolverConfig *cfg, double *root) {
    double bound = root_bound(p);
    const int pairs = 12;

    for (int i = 0; i < pairs; ++i) {
        double x0 = -bound + (2.0 * bound * i) / pairs;
        double x1 = x0 + bound / 7.0;

        for (int it = 0; it < cfg->max_iter; ++it) {
            double f0 = eval_poly(p, x0);
            double f1 = eval_poly(p, x1);
            double denom = f1 - f0;
            if (fabs(denom) < 1e-14) {
                break;
            }
            double x2 = x1 - f1 * (x1 - x0) / denom;
            if (!isfinite(x2)) {
                break;
            }
            if (fabs(x2 - x1) < cfg->tol_x && fabs(eval_poly(p, x2)) < cfg->tol_f * 100.0) {
                *root = x2;
                return true;
            }
            x0 = x1;
            x1 = x2;
        }
    }
    return false;
}

static bool sample_refine_method(const Polynomial *p, const SolverConfig *cfg, double *root) {
    double bound = root_bound(p);
    int steps = cfg->scan_steps > 500 ? cfg->scan_steps : 500;
    double best_x = -bound;
    double best_abs_f = fabs(eval_poly(p, best_x));

    for (int i = 1; i <= steps; ++i) {
        double x = -bound + (2.0 * bound * i) / steps;
        double f = fabs(eval_poly(p, x));
        if (f < best_abs_f) {
            best_abs_f = f;
            best_x = x;
        }
    }

    if (best_abs_f > 1e-3) {
        return false;
    }

    double x = best_x;
    for (int i = 0; i < cfg->max_iter; ++i) {
        double fx = eval_poly(p, x);
        double dfx = eval_poly_derivative(p, x);
        if (fabs(dfx) < 1e-14) {
            break;
        }
        double nx = x - fx / dfx;
        if (!isfinite(nx)) {
            break;
        }
        if (fabs(nx - x) < cfg->tol_x && fabs(eval_poly(p, nx)) < cfg->tol_f * 100.0) {
            *root = nx;
            return true;
        }
        x = nx;
    }

    return false;
}

static void horner_deflate(Polynomial *p, double root) {
    int new_degree = p->degree - 1;
    double next[MAX_DEGREE + 1] = {0};
    next[0] = p->coeffs[0];
    for (int i = 1; i <= new_degree; ++i) {
        next[i] = p->coeffs[i] + next[i - 1] * root;
    }

    for (int i = 0; i <= new_degree; ++i) {
        p->coeffs[i] = next[i];
    }
    p->degree = new_degree;
    remove_leading_zeros(p, 1e-10);
}

static void deduplicate_and_round(RealRoots *roots) {
    const double dup_eps = 1e-4;
    const double int_eps = 1e-2;

    for (size_t i = 0; i < roots->count; ++i) {
        for (size_t j = i + 1; j < roots->count; ++j) {
            if (roots->values[j] < roots->values[i]) {
                double t = roots->values[i];
                roots->values[i] = roots->values[j];
                roots->values[j] = t;
            }
        }
    }

    double unique[6] = {0};
    size_t n = 0;
    for (size_t i = 0; i < roots->count; ++i) {
        if (n == 0 || fabs(roots->values[i] - unique[n - 1]) > dup_eps) {
            unique[n++] = roots->values[i];
        }
    }

    roots->count = n;
    for (size_t i = 0; i < n; ++i) {
        double v = unique[i];
        double iv = nearbyint(v);
        if (fabs(v - iv) <= int_eps) {
            roots->values[i] = iv;
        } else {
            roots->values[i] = round(v * 1e6) / 1e6;
        }
    }
}

SolveStatus solve_equation(const SexticCoefficients *coeffs,
                           const SolverConfig *cfg,
                           RealRoots *out_roots,
                           char *err_buf,
                           size_t err_size) {
    if (!validate_coefficients(coeffs, err_buf, err_size)) {
        return SOLVE_ERROR;
    }
    if (out_roots == NULL || cfg == NULL) {
        snprintf(err_buf, err_size, "Internal error: null output/config pointer.");
        return SOLVE_ERROR;
    }

    out_roots->count = 0;

    Polynomial p = {.coeffs = {coeffs->a, coeffs->b, coeffs->c, coeffs->d, coeffs->n, coeffs->m, coeffs->k}, .degree = 6};

    bool all_zero = true;
    for (int i = 0; i <= 6; ++i) {
        if (fabs(p.coeffs[i]) > EPS_ZERO) {
            all_zero = false;
            break;
        }
    }
    if (all_zero) {
        return SOLVE_ANY_NUMBER;
    }

    remove_leading_zeros(&p, EPS_ZERO);

    if (p.degree == 0) {
        if (fabs(p.coeffs[0]) <= EPS_ZERO) {
            return SOLVE_ANY_NUMBER;
        }
        return SOLVE_NO_SOLUTIONS;
    }

    while (p.degree >= 4) {
        double root = 0.0;
        bool found = newton_method(&p, cfg, &root);
        if (!found) {
            found = chord_tangent_method(&p, cfg, &root);
        }
        if (!found) {
            found = secant_method(&p, cfg, &root);
        }
        if (!found) {
            found = sample_refine_method(&p, cfg, &root);
        }

        if (!found || !isfinite(root)) {
            if ((p.degree % 2) == 0) {
                break;
            }
            snprintf(
                err_buf,
                err_size,
                "Failed to find a real root for odd degree %d polynomial using Newton/chord-tangent/secant methods.",
                p.degree);
            return SOLVE_ERROR;
        }

        if (!add_root(out_roots, root)) {
            snprintf(err_buf, err_size, "Internal error: too many roots collected.");
            return SOLVE_ERROR;
        }

        horner_deflate(&p, root);
    }

    if (p.degree >= 4) {
        if (out_roots->count == 0) {
            return SOLVE_NO_SOLUTIONS;
        }
        deduplicate_and_round(out_roots);
        return SOLVE_OK;
    }

    if (p.degree == 3) {
        double cubic_roots[3] = {0};
        size_t rc = cardano_method(&p, cubic_roots);
        for (size_t i = 0; i < rc; ++i) {
            if (!add_root(out_roots, cubic_roots[i])) {
                snprintf(err_buf, err_size, "Internal error: too many roots collected.");
                return SOLVE_ERROR;
            }
        }
    } else if (p.degree == 2) {
        double qr[2] = {0};
        size_t rc = vieta_method(&p, qr);
        for (size_t i = 0; i < rc; ++i) {
            if (!add_root(out_roots, qr[i])) {
                snprintf(err_buf, err_size, "Internal error: too many roots collected.");
                return SOLVE_ERROR;
            }
        }
    } else if (p.degree == 1) {
        double lr = 0.0;
        if (linear_method(&p, &lr)) {
            if (!add_root(out_roots, lr)) {
                snprintf(err_buf, err_size, "Internal error: too many roots collected.");
                return SOLVE_ERROR;
            }
        }
    } else {
        if (fabs(p.coeffs[0]) <= EPS_ZERO) {
            return SOLVE_ANY_NUMBER;
        }
        return SOLVE_NO_SOLUTIONS;
    }

    deduplicate_and_round(out_roots);
    if (out_roots->count == 0) {
        return SOLVE_NO_SOLUTIONS;
    }
    return SOLVE_OK;
}
