#!/bin/bash

# All the constant options here:

OPTIONS="--init-motion models/satellite3D_init.rkx"
OPTIONS="${OPTIONS} --inertia models/sat3D_airship_inertia.rkx"
OPTIONS="${OPTIONS} --Q-matrix models/sat3D_airship_Q.rkx"
OPTIONS="${OPTIONS} --IMU-config models/sat3D_airship_IMU_config.rkx"
OPTIONS="${OPTIONS} --time-step 0.01"
# OPTIONS="${OPTIONS} --imkf"
# OPTIONS="${OPTIONS} --imkf-em"
OPTIONS="${OPTIONS} --imkf-emd"
# OPTIONS="${OPTIONS} --imkf-emdJ"

# OPTIONS="${OPTIONS} --R-matrix temp_sat3D_airship_R.rkx"

OPTIONS="${OPTIONS} --R-matrix temp_sat3D_airship_R_gyro.rkx"
OPTIONS="${OPTIONS} --gyro"

OPTIONS="${OPTIONS} --prediction-runs --prediction-interval 5.0"

OPTIONS="${OPTIONS} --pred-assumption 0"
# OPTIONS="${OPTIONS} --pred-assumption 1"

# OPTIONS="${OPTIONS} --tsosakf --Pa-matrix temp_sat3D_airship_Pa_em.rkx"
OPTIONS="${OPTIONS} --tsosakf --Pa-matrix temp_sat3D_airship_Pa_emd.rkx"
# OPTIONS="${OPTIONS} --tsosakf --Pa-matrix temp_sat3D_airship_Pa_emdJ.rkx"

OPTIONS="${OPTIONS} --monte-carlo --min-skips 4 --max-skips 6"


DM_COV=1e-6
ECC_COV=1e-6
LD_COV=20
RD_COV=1

# POS_MEAS_COV=0.0005
# ANG_MEAS_COV=0.0009
# GYRO_MEAS_COV=0.0005
POS_MEAS_COV=0.00005
ANG_MEAS_COV=0.00009
GYRO_MEAS_COV=0.00005

echo "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\" ?>
<!DOCTYPE reak_serialization>
<reak_serialization version=\"2\">
<steady_param_covariance type_ID=\"18.4.5.1.0\" version=\"1\">
    <q_count>\"6\"</q_count>
    <q_q[0]>\"$DM_COV\"</q_q[0]>
    <q_q[1]>\"$ECC_COV\"</q_q[1]>
    <q_q[2]>\"$ECC_COV\"</q_q[2]>
    <q_q[3]>\"$ECC_COV\"</q_q[3]>
    <q_q[4]>\"$LD_COV\"</q_q[4]>
    <q_q[5]>\"$RD_COV\"</q_q[5]>
    <rowCount>\"6\"</rowCount>
</steady_param_covariance>
</reak_serialization>
" > temp_sat3D_airship_Pa_emd.rkx

echo "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\" ?>
<!DOCTYPE reak_serialization>
<reak_serialization version=\"2\">
<measurement_noise type_ID=\"18.4.5.1.0\" version=\"1\">
    <q_count>\"9\"</q_count>
    <q_q[0]>\"$POS_MEAS_COV\"</q_q[0]>
    <q_q[1]>\"$POS_MEAS_COV\"</q_q[1]>
    <q_q[2]>\"$POS_MEAS_COV\"</q_q[2]>
    <q_q[3]>\"$ANG_MEAS_COV\"</q_q[3]>
    <q_q[4]>\"$ANG_MEAS_COV\"</q_q[4]>
    <q_q[5]>\"$ANG_MEAS_COV\"</q_q[5]>
    <q_q[6]>\"$GYRO_MEAS_COV\"</q_q[6]>
    <q_q[7]>\"$GYRO_MEAS_COV\"</q_q[7]>
    <q_q[8]>\"$GYRO_MEAS_COV\"</q_q[8]>
    <rowCount>\"9\"</rowCount>
</measurement_noise>
</reak_serialization>
" > temp_sat3D_airship_R_gyro.rkx

IN_OPTIONS="--input $1 --input-format ssv"
OUT_OPTIONS="--output $1 --output-format ssv"
./estimate_satellite3D ${OPTIONS} ${IN_OPTIONS} ${OUT_OPTIONS}




