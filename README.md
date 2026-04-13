# Sextic Equation Solver (C)

CLI program for equations of the form:

```
a x^6 + b x^5 + c x^4 + d x^3 + n x^2 + m x + k = 0
```

## Features

- Secure interactive input: requires exactly 7 numeric coefficients, otherwise re-prompts.
- Handles special cases:
  - all coefficients are `0` -> `Any number is a solution.`
  - reduced constant equation `c = 0` -> `Any number.`
  - reduced constant equation `c != 0` -> `No solutions.`
- Removes leading zero coefficients automatically (degree reduction).
- Solves by degree:
  - degree 1: `linear_method`
  - degree 2: `vieta_method`
  - degree 3: `cardano_method`
  - degree >= 4: numerical methods in priority order:
    1) `newton_method`
    2) `chord_tangent_method`
    3) `secant_method`
- For degree >= 4, applies Horner deflation after each found root.
- Post-processing:
  - duplicate root removal (`1e-4`)
  - integer snapping (`1e-2`)
  - rounding to 6 decimals
  - output format: `Roots: [x1, x2, ...]`
- Repeat loop with `Continue? (y/n)`.

## Build

```bash
make
```

## Run

```bash
./sextic_solver
```

## Tests

```bash
make test           # quick smoke test
make test-diverse   # normal, edge, and numerical cases
```
