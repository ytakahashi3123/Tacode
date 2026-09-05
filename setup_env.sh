#!/bin/bash
#
# Tacode -- optional environment setup
#
# Builds a Python virtual environment with the packages Tacode needs, and
# verifies that they import. Running Tacode does NOT require this script: if you
# already have numpy, scipy, PyYAML and simplekml available, keep using your own
# interpreter and ignore this.
#
# Once the environment exists, run_tacode.sh picks it up automatically.
#
#   ./setup_env.sh                 build (or reuse) ./.venv and install
#   ./setup_env.sh --check         only report what is missing, install nothing
#   ./setup_env.sh --upgrade       upgrade the packages in an existing venv
#   ./setup_env.sh --with-optional also install gpxpy
#   ./setup_env.sh --venv <path>   use a different location for the venv
#   ./setup_env.sh --python <cmd>  base interpreter (default: python3)
#   ./setup_env.sh --help
#

set -u

TACODE_ROOT=$(cd "$(dirname "$0")" && pwd)
REQ_FILE=$TACODE_ROOT/requirements.txt

VENV_DIR=$TACODE_ROOT/.venv
BASE_PYTHON=python3
MODE=install
WITH_OPTIONAL=0

REQUIRED_MODULES="numpy scipy yaml simplekml"
OPTIONAL_MODULES="gpxpy"

usage() {
  # ファイル先頭のコメントブロック（shebang の次から最初の非コメント行まで）をそのまま出す
  awk 'NR == 1 { next } /^#/ { sub(/^# ?/, ""); print; next } { exit }' "$0"
}

while [ $# -gt 0 ]; do
  case "$1" in
    --check)         MODE=check ;;
    --upgrade)       MODE=upgrade ;;
    --with-optional) WITH_OPTIONAL=1 ;;
    --venv)          shift; [ $# -gt 0 ] || { echo "--venv needs a path" >&2; exit 2; }; VENV_DIR=$1 ;;
    --python)        shift; [ $# -gt 0 ] || { echo "--python needs a command" >&2; exit 2; }; BASE_PYTHON=$1 ;;
    -h|--help)       usage; exit 0 ;;
    *)               echo "Unknown option: $1" >&2; echo "Try: $0 --help" >&2; exit 2 ;;
  esac
  shift
done

# ---------------------------------------------------------------------------
# report_modules <python> <modules...>
#   Prints one line per module and returns the number that are missing.
# ---------------------------------------------------------------------------
report_modules() {
  local py=$1; shift
  local missing=0
  local m
  for m in "$@"; do
    local ver
    if ver=$("$py" -c "import $m,sys; print(getattr($m,'__version__','(no __version__)'))" 2>/dev/null); then
      printf '  %-12s %s\n' "$m" "$ver"
    else
      printf '  %-12s MISSING\n' "$m"
      missing=$((missing + 1))
    fi
  done
  return $missing
}

# ---------------------------------------------------------------------------
# --check : inspect without changing anything
# ---------------------------------------------------------------------------
if [ "$MODE" = check ]; then
  if [ -x "$VENV_DIR/bin/python" ]; then
    CHECK_PYTHON=$VENV_DIR/bin/python
    echo "Checking the Tacode virtual environment: $VENV_DIR"
  elif command -v "$BASE_PYTHON" > /dev/null 2>&1; then
    CHECK_PYTHON=$BASE_PYTHON
    echo "No virtual environment at $VENV_DIR"
    echo "Checking $BASE_PYTHON instead: $(command -v "$BASE_PYTHON")"
  else
    echo "Neither $VENV_DIR nor '$BASE_PYTHON' is available." >&2
    exit 1
  fi

  echo "--$("$CHECK_PYTHON" --version 2>&1)"
  echo "Required:"
  report_modules "$CHECK_PYTHON" $REQUIRED_MODULES
  missing=$?
  echo "Optional:"
  report_modules "$CHECK_PYTHON" $OPTIONAL_MODULES

  if [ "$missing" -gt 0 ]; then
    echo
    echo "$missing required package(s) missing. Run: $0"
    exit 1
  fi
  echo
  echo "All required packages are available."
  exit 0
fi

# ---------------------------------------------------------------------------
# install / upgrade
# ---------------------------------------------------------------------------
if [ ! -f "$REQ_FILE" ]; then
  echo "Requirements file not found: $REQ_FILE" >&2
  exit 1
fi

if ! command -v "$BASE_PYTHON" > /dev/null 2>&1; then
  echo "Base interpreter not found: $BASE_PYTHON" >&2
  echo "Specify one with: $0 --python /path/to/python3" >&2
  exit 1
fi

if [ -x "$VENV_DIR/bin/python" ]; then
  echo "Reusing the existing virtual environment: $VENV_DIR"
else
  echo "Creating a virtual environment: $VENV_DIR"
  if ! "$BASE_PYTHON" -m venv "$VENV_DIR"; then
    echo >&2
    echo "Failed to create the virtual environment." >&2
    echo "On Debian/Ubuntu the venv module ships separately:" >&2
    echo "  sudo apt install python3-venv" >&2
    exit 1
  fi
fi

VENV_PYTHON=$VENV_DIR/bin/python
echo "--$("$VENV_PYTHON" --version 2>&1)"

PIP_FLAGS=""
if [ "$MODE" = upgrade ]; then
  PIP_FLAGS="--upgrade"
fi

echo "Installing the required packages..."
if ! "$VENV_PYTHON" -m pip install --quiet $PIP_FLAGS -r "$REQ_FILE"; then
  echo "Installation failed. Re-run without --quiet to see the details:" >&2
  echo "  $VENV_PYTHON -m pip install -r $REQ_FILE" >&2
  exit 1
fi

if [ "$WITH_OPTIONAL" -eq 1 ]; then
  echo "Installing the optional packages..."
  if ! "$VENV_PYTHON" -m pip install --quiet $PIP_FLAGS $OPTIONAL_MODULES; then
    echo "Failed to install the optional packages: $OPTIONAL_MODULES" >&2
    exit 1
  fi
fi

echo "Verifying..."
echo "Required:"
report_modules "$VENV_PYTHON" $REQUIRED_MODULES
missing=$?
if [ "$missing" -gt 0 ]; then
  echo >&2
  echo "$missing required package(s) still missing after installation." >&2
  exit 1
fi

if [ "$WITH_OPTIONAL" -eq 1 ]; then
  echo "Optional:"
  report_modules "$VENV_PYTHON" $OPTIONAL_MODULES
fi

cat <<MSG

Done. The environment is ready:
  $VENV_PYTHON

run_tacode.sh now picks it up automatically, so you can just run:
  cd tutorial/work
  ./run_tacode.sh

To use a different interpreter for one run, set TACODE_PYTHON:
  TACODE_PYTHON=/path/to/python ./run_tacode.sh
MSG
