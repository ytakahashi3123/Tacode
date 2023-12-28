#!/bin/bash

TACODE_HOME=$1
LD=$TACODE_HOME/tacode.py
LOG=log_tacode

python3.9 $LD > $LOG
#python3.9 $LD 2>&1 | tee $LOG
