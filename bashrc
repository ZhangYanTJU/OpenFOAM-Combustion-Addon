#!/bin/bash

FOAM_OFCA_DIR=$WM_PROJECT_INST_DIR/CombustionAddons
alias ofca="cd $FOAM_OFCA_DIR"

#CUBA_PATH=$HOME/Programme/Cuba-3.0/platforms/$WM_OPTIONS
CUBA_PATH=/usr
export CUBA_INC_PATH=$CUBA_PATH/include
export CUBA_LIB_PATH=$CUBA_PATH/lib

CANTERA_PATH=/usr/local/cantera
export CANTERA_INC_PATH=$CANTERA_PATH/include
export CANTERA_LIB_PATH=$CANTERA_PATH/lib

export OCTAVE_INC_PATH=/usr/include/octave-3.2.4
export OCTAVE_LIB_PATH=/usr/lib/octave-3.2.4

cleanProg=$WM_PROJECT_DIR/bin/foamCleanPath

if [ -d $FOAM_OFCA_DIR/bin/$WM_OPTIONS ]; then
 export PATH=$FOAM_OFCA_DIR/bin/$WM_OPTIONS:$PATH 
 cleanEnv=`$cleanProg "$PATH"` && PATH="$cleanEnv"
fi
if [ -d $FOAM_OFCA_DIR/lib/$WM_OPTIONS ]; then
 export LD_LIBRARY_PATH=$OCTAVE_LIB_PATH:$FOAM_OFCA_DIR/lib/$WM_OPTIONS:$LD_LIBRARY_PATH
 cleanEnv=`$cleanProg "$LD_LIBRARY_PATH"` && LD_LIBRARY_PATH="$cleanEnv"
fi

