#!/bin/bash

BASE_MDL_NAME="models/satellite3D"

MAIN_OPTIONS="--init=\"${BASE_MDL_NAME}_init.rkx\""
MAIN_OPTIONS="${MAIN_OPTIONS} --inertia=\"${BASE_MDL_NAME}_inertia.rkx\""
MAIN_OPTIONS="${MAIN_OPTIONS} --Q-matrix=\"${BASE_MDL_NAME}_Q.rkx\""
MAIN_OPTIONS="${MAIN_OPTIONS} --R-matrix=\"${BASE_MDL_NAME}_R.rkx\""
MAIN_OPTIONS="${MAIN_OPTIONS} --R-added=\"${BASE_MDL_NAME}_R_added.rkx\""
MAIN_OPTIONS="${MAIN_OPTIONS} --IMU-config=\"${BASE_MDL_NAME}_IMU_config.rkx\""

MAIN_OPTIONS="${MAIN_OPTIONS} --system-output=\"${BASE_MDL_NAME}\""

EXTRA_OPTIONS="--time-step=0.01"

echo "./estimate_satellite3D --generate-mdl-files ${MAIN_OPTIONS} ${EXTRA_OPTIONS}"
./estimate_satellite3D --generate-mdl-files ${MAIN_OPTIONS} ${EXTRA_OPTIONS}
