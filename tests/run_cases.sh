#!/usr/bin/env bash
set -euo pipefail

BIN="${1:-./sextic_solver}"

run_case() {
  local category="$1"
  local name="$2"
  local input="$3"
  local expected="$4"

  printf '=== [%s] %s ===\n' "$category" "$name"
  local out
  out=$(printf "%b" "$input" | "$BIN")
  printf '%s\n' "$out"

  if [[ "$out" == *"$expected"* ]]; then
    printf 'PASS\n\n'
  else
    printf 'FAIL: expected to find "%s"\n\n' "$expected" >&2
    return 1
  fi
}

# Normal cases with exact roots
run_case "normal" "six exact integer roots" "1 -21 175 -735 1624 -1764 720\nn\n" "Roots: [1, 2, 3, 4, 5, 6]"
run_case "normal" "two exact integer roots" "1 0 0 0 0 0 -1\nn\n" "Roots: [-1, 1]"
run_case "normal" "cubic exact roots" "0 0 0 1 -6 11 -6\nn\n" "Roots: [1, 2, 3]"

# Edge cases
run_case "edge" "all coefficients are zero" "0 0 0 0 0 0 0\nn\n" "Any number is a solution."
run_case "edge" "constant-only nonzero" "0 0 0 0 0 0 5\nn\n" "No solutions."
run_case "edge" "invalid input then valid" "x y z\n1 0 0 0 0 0 -1\nn\n" "Invalid input. Please enter exactly 7 numeric values."

# Numerical cases
run_case "numerical" "repeated real roots and complex pair" "1 2 -2 -2 1 -4 4\nn\n" "Roots: [-2, 1]"
run_case "numerical" "no real roots sextic" "1 0 0 0 0 0 1\nn\n" "No solutions."
run_case "numerical" "non-integer exact roots" "1 -5.05 -5.95 46.4815 -30.4695 -16.56 5.346\nn\n" "Roots: [-3, -0.500000, 0.250000, 1.200000, 2.700000, 4.400000]"

printf 'All requested test categories passed.\n'
