#!/bin/bash

OPTIONS="--init-motion models/satellite3D_init.rkx"
OPTIONS="${OPTIONS} --inertia models/sat3D_airship_inertia.rkx"
OPTIONS="${OPTIONS} --Q-matrix models/sat3D_airship_Q.rkx"
OPTIONS="${OPTIONS} --IMU-config models/sat3D_airship_IMU_config.rkx"
OPTIONS="${OPTIONS} --time-step 0.01"
# OPTIONS="${OPTIONS} --imkf"
# OPTIONS="${OPTIONS} --imkf-em"
OPTIONS="${OPTIONS} --imkf-emd"

# OPTIONS="${OPTIONS} --R-matrix models/sat3D_airship_R.rkx"

OPTIONS="${OPTIONS} --R-matrix models/sat3D_airship_R_gyro.rkx"
OPTIONS="${OPTIONS} --gyro"

OPTIONS="${OPTIONS} --prediction-runs --prediction-interval 5.0"

OPTIONS="${OPTIONS} --pred-assumption 0"
# OPTIONS="${OPTIONS} --pred-assumption 1"


# IN_OPTIONS="--input $1 --input-format tsv"
# OUT_OPTIONS="--output $1 --output-format ssv"
# ./estimate_satellite3D ${OPTIONS} ${IN_OPTIONS} ${OUT_OPTIONS}


for f in $1/*.tsv
do
  IN_OPTIONS="--input $f --input-format tsv"
  OUT_OPTIONS="--output $f --output-format ssv"
  ./estimate_satellite3D ${OPTIONS} ${IN_OPTIONS} ${OUT_OPTIONS}
done
