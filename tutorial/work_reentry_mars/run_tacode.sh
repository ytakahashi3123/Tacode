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

# The Mars atmosphere table is not in the repository: it is Mars Climate Database data,
# which each user fetches under the MCD terms. Stop with the command if it is missing.
TABLE=database/atmosphere/atmospheremodel_mars.txt
if [ ! -f "$TABLE" ]; then
  echo "The Mars atmosphere table is missing: $TABLE"
  echo "It comes from the Mars Climate Database (MCD, LMD/OU/IAA/ESA/CNES) and is not"
  echo "distributed with Tacode. Fetch it from the MCD web interface (network required):"
  echo ""
  echo "  $PYTHON_RUN ../../database/atmosphere/generate_atmosphere_table_mars.py --latitude 19.13 \\"
  echo "      --longitude -33.22 --ls 142.7 --local-time 3.0 --altitude 0 200 \\"
  echo "      -o $TABLE"
  echo ""
  echo "MCD terms (https://www-mars.lmd.jussieu.fr/mars/access.html): scientific use is free"
  echo "provided the origin of the data is quoted and the MCD team is kept informed;"
  echo "no commercial use without their authorization."
  exit 1
fi
LOG=log_tacode

export OMP_NUM_THREADS=1

$PYTHON_RUN $LD > $LOG
#$PYTHON_RUN  $LD 2>&1 | tee $LOG