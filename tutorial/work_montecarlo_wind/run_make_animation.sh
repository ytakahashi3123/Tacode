#!/bin/bash

TACODE_HELPER_HOME=../../src_helper

if [ -n "${TACODE_PYTHON:-}" ]; then
  PYTHON_RUN=$TACODE_PYTHON
elif [ -x "$TACODE_HOME/../.venv/bin/python" ]; then
  PYTHON_RUN=$TACODE_HOME/../.venv/bin/python
else
  PYTHON_RUN=python3
fi

LD=$TACODE_HELPER_HOME/montecarlo_animation/montecarlo_animation.py
LOG=log_montecarlo_animation

$PYTHON_RUN $LD work_montecarlo_wind 2>&1 | tee $LOG