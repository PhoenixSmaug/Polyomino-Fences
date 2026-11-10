#!/bin/bash

# ── Configuration ────────────────────────────────────────────────────────────
TIMEOUT=60          # seconds before we consider the run "too slow"
LOWER=19000             # inclusive lower bound for the search
UPPER=20000         # inclusive upper bound for the search
FIXED_ARG=7         # n_max
# ─────────────────────────────────────────────────────────────────────────────

SCRIPT="best-fence-estimate-hexa-step1.py"

echo "Binary search: smallest value in [$LOWER, $UPPER] that finishes within ${TIMEOUT}s"
echo "Command: timeout ${TIMEOUT}s pypy3 $SCRIPT $FIXED_ARG <value>"
echo ""

# Quick sanity checks
echo "Checking bounds..."

timeout "${TIMEOUT}s" pypy3 "$SCRIPT" "$FIXED_ARG" "$UPPER"
if [ $? -ne 0 ]; then
    echo "ERROR: Upper bound $UPPER does not finish within ${TIMEOUT}s. Raise UPPER or TIMEOUT."
    exit 1
fi
echo "  Upper bound $UPPER: OK (finishes in time)"

timeout "${TIMEOUT}s" pypy3 "$SCRIPT" "$FIXED_ARG" "$LOWER"
if [ $? -eq 0 ]; then
    echo "  Lower bound $LOWER: also finishes in time — smallest possible value is $LOWER (or lower)."
    echo "Result: $LOWER"
    exit 0
fi
echo "  Lower bound $LOWER: too slow (expected)"
echo ""

# Binary search
#   Invariant: LOWER is always "too slow", UPPER is always "fast enough"
lo=$LOWER
hi=$UPPER
best=$UPPER

while [ $((hi - lo)) -gt 1 ]; do
    mid=$(( (lo + hi) / 2 ))
    echo -n "  Testing value=$mid ... "

    timeout "${TIMEOUT}s" pypy3 "$SCRIPT" "$FIXED_ARG" "$mid"
    status=$?

    if [ $status -eq 0 ]; then
        echo "OK  → search lower half [$lo, $mid]"
        best=$mid
        hi=$mid
    else
        echo "SLOW (exit $status) → search upper half [$mid, $hi]"
        lo=$mid
    fi
done

echo ""
echo "══════════════════════════════════════════════"
echo "Smallest value that finishes within ${TIMEOUT}s: $best"
echo "══════════════════════════════════════════════"
