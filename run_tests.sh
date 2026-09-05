#!/bin/bash
#
# Tacode -- run the test suite
#
# The tests use only the standard library's unittest, so they need nothing
# beyond the packages Tacode itself requires (numpy, scipy, PyYAML).
#
#   ./run_tests.sh                     run everything
#   ./run_tests.sh -v                  verbose
#   ./run_tests.sh test_kepler         run one module
#   ./run_tests.sh test_kepler -v      both
#
# The interpreter is resolved the same way run_tacode.sh does:
#   1. $TACODE_PYTHON  2. ./.venv/bin/python  3. python3
#

set -u

TACODE_ROOT=$(cd "$(dirname "$0")" && pwd)

if [ -n "${TACODE_PYTHON:-}" ]; then
  PYTHON_RUN=$TACODE_PYTHON
elif [ -x "$TACODE_ROOT/.venv/bin/python" ]; then
  PYTHON_RUN=$TACODE_ROOT/.venv/bin/python
else
  PYTHON_RUN=python3
fi

VERBOSE=""
MODULES=""
for arg in "$@"; do
  case "$arg" in
    -v|--verbose) VERBOSE="-v" ;;
    -h|--help)    awk 'NR == 1 { next } /^#/ { sub(/^# ?/, ""); print; next } { exit }' "$0"; exit 0 ;;
    *)            MODULES="$MODULES $arg" ;;
  esac
done

cd "$TACODE_ROOT/test" || exit 1

if [ -n "$MODULES" ]; then
  # shellcheck disable=SC2086
  exec "$PYTHON_RUN" -m unittest $VERBOSE $MODULES
fi

# shellcheck disable=SC2086
exec "$PYTHON_RUN" -m unittest discover -s . -p 'test_*.py' $VERBOSE
