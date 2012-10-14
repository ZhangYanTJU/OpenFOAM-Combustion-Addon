#!/bin/bash

OFCADIR=$WM_PROJECT_INST_DIR/CombustionAddons
alias ofca="cd $OFCADIR"

#CUBA_PATH=$HOME/Programme/Cuba-3.0/platforms/$WM_OPTIONS
CUBA_PATH=/usr
export CUBA_INC_PATH=$CUBA_PATH/include
export CUBA_LIB_PATH=$CUBA_PATH/lib

CANTERA_PATH=/usr/local/cantera
export CANTERA_INC_PATH=$CANTERA_PATH/include
export CANTERA_LIB_PATH=$CANTERA_PATH/lib
