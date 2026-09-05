#!/bin/bash

TACODE_HOME=../../src

# Python interpreter, in order of preference:
#   1. $TACODE_PYTHON, if set
#   2. the virtual environment built by setup_env.sh, if it exists
#   3. python3
if [ -n "${TACODE_PYTHON:-}" ]; then
  PYTHON_RUN=$TACODE_PYTHON
elif [ -x "$TACODE_HOME/../.venv/bin/python" ]; then
  PYTHON_RUN=$TACODE_HOME/../.venv/bin/python
else
  PYTHON_RUN=python3
fi
LD=$TACODE_HOME/tacode.py
LOG=log_tacode

export OMP_NUM_THREADS=1

$PYTHON_RUN $LD > $LOG
#$PYTHON_RUN  $LD 2>&1 | tee $LOG