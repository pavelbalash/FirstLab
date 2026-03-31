# Sextic Equation Solver (C Skeleton)

Project skeleton for a robust CLI program that solves sixth-degree equations of the form:

```
a x^6 + b x^5 + c x^4 + d x^3 + n x^2 + m x + k = 0
```

## Current status

This is a **starter skeleton** with:
- Input parsing
- Coefficient validation
- Horner polynomial evaluation
- Derivative evaluation
- Cauchy root bound utility
- Placeholder `find_real_roots()` for future implementation

## Build

```bash
make
```

## Run

```bash
./sextic_solver
```

## Demo mode

```bash
./sextic_solver --demo
```

## Next implementation tasks

1. Interval isolation in `[-R, R]` where `R` is Cauchy bound.
2. Bisection root refinement.
3. Optional Newton acceleration with derivative safety checks.
4. Duplicate root suppression with tolerance.
5. Comprehensive tests.
