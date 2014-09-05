#!/bin/bash


#####################################################################
# Evaluate a floating point number expression.

function float_eval()
{
    local stat=0
    local result=0.0
    if [[ $# -gt 0 ]]; then
        result=$(echo "scale=20; $*" | bc -q 2>/dev/null)
        stat=$?
        if [[ $stat -eq 0  &&  -z "$result" ]]; then stat=1; fi
    fi
    echo $result
    return $stat
}


#####################################################################
# Evaluate a floating point number conditional expression.

function float_cond()
{
    local cond=0
    if [[ $# -gt 0 ]]; then
        cond=$(echo "$*" | bc -q 2>/dev/null)
        if [[ -z "$cond" ]]; then cond=0; fi
        if [[ "$cond" != 0  &&  "$cond" != 1 ]]; then cond=0; fi
    fi
    local stat=$((cond == 0))
    return $stat
}


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

if [[ "0" == "0" ]]
then

INTERVAL_DIV=4.0

DM_COV_LO=0.0000001
DM_COV_HI=0.0001

ECC_COV_LO=0.0000001
ECC_COV_HI=0.0001

LD_COV_LO=0.1
LD_COV_HI=100

RD_COV_LO=0.1
RD_COV_HI=100

POS_MEAS_COV_LO=0.000005
POS_MEAS_COV_HI=0.0005
ANG_MEAS_COV_LO=0.000009
ANG_MEAS_COV_HI=0.0009
GYRO_MEAS_COV_LO=0.000005
GYRO_MEAS_COV_HI=0.0005

DM_COV=$DM_COV_HI
ECC_COV=$ECC_COV_HI
LD_COV=$LD_COV_HI

POS_MEAS_COV=$POS_MEAS_COV_HI
ANG_MEAS_COV=$ANG_MEAS_COV_HI
GYRO_MEAS_COV=$GYRO_MEAS_COV_HI

FINAL_OUTPUT_FILE="$1/all_predstats.ssv"
echo "% pos_meas_cov ang_meas_cov gyro_meas_cov dm_cov ecc_cov ld_cov rd_cov P_th success_rate pred_start_time ep_m ea_m ev_m ew_m pdf_est lr_est pdf_meas lr_meas" > $FINAL_OUTPUT_FILE

(( MC_TRIALS_COUNT = 0 ))

POS_MEAS_COV=$POS_MEAS_COV_HI
until (( $(echo "$POS_MEAS_COV > $POS_MEAS_COV_LO" | bc -q 2>/dev/null) == "0" ))
do
  echo "pos-meas cov = $POS_MEAS_COV"
  ANG_MEAS_COV=$ANG_MEAS_COV_HI
  until (( $(echo "$ANG_MEAS_COV > $ANG_MEAS_COV_LO" | bc -q 2>/dev/null) == "0" ))
  do
    echo "ang-meas cov = $ANG_MEAS_COV"
    GYRO_MEAS_COV=$GYRO_MEAS_COV_HI
    until (( $(echo "$GYRO_MEAS_COV > $GYRO_MEAS_COV_LO" | bc -q 2>/dev/null) == "0" ))
    do
      echo "gyro-meas cov = $GYRO_MEAS_COV"
      DM_COV=$DM_COV_HI
      until (( $(echo "$DM_COV > $DM_COV_LO" | bc -q 2>/dev/null) == "0" ))
      do
        echo "dm cov = $DM_COV"
        ECC_COV=$ECC_COV_HI
        until (( $(echo "$ECC_COV > $ECC_COV_LO" | bc -q 2>/dev/null) == "0" ))
        do
          echo "ecc. cov = $ECC_COV"
          LD_COV=$LD_COV_HI
          until (( $(echo "$LD_COV > $LD_COV_LO" | bc -q 2>/dev/null) == "0" ))
          do
            echo "lin-drag cov = $LD_COV"
            RD_COV=$RD_COV_HI
            until (( $(echo "$RD_COV > $RD_COV_LO" | bc -q 2>/dev/null) == "0" ))
            do
              echo "rot-drag cov = $RD_COV"
              
              if [[ "0" == "1" ]]
              then
              ACCUM_FILE="$1/predstats_${POS_MEAS_COV}_${ANG_MEAS_COV}_${GYRO_MEAS_COV}_${DM_COV}_${ECC_COV}_${LD_COV}_${RD_COV}.ssv"
              echo "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0" > $ACCUM_FILE
              
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
              
              for f in $1/*
              do
                IN_OPTIONS="--input \"$f\" --input-format ssv"
                OUT_OPTIONS="--output $1/tmp_predstats.ssv --output-format ssv"
                ./estimate_satellite3D ${OPTIONS} ${IN_OPTIONS} ${OUT_OPTIONS}
                for g in $1/tmp_predstats*.ssv
                do
                  octave --quiet --norc --no-window-system --eval "add_predstats('$g','$ACCUM_FILE')"
                done
              done
              
              # Put all accumulated data into a single file:
              while read p; do
                echo "POS_MEAS_COV ANG_MEAS_COV GYRO_MEAS_COV DM_COV ECC_COV LD_COV RD_COV $p" >> $FINAL_OUTPUT_FILE
              done <$ACCUM_FILE
              
              else
              
              (( MC_TRIALS_COUNT++ ))
              
              fi
              
              RD_COV=$(echo "scale=20; $RD_COV / $INTERVAL_DIV" | bc -q 2>/dev/null)
            done
            LD_COV=$(float_eval "$LD_COV / $INTERVAL_DIV")
          done
          ECC_COV=$(float_eval "$ECC_COV / $INTERVAL_DIV")
        done
        DM_COV=$(float_eval "$DM_COV / $INTERVAL_DIV")
      done
      GYRO_MEAS_COV=$(float_eval "$GYRO_MEAS_COV / $INTERVAL_DIV")
    done
    ANG_MEAS_COV=$(float_eval "$ANG_MEAS_COV / $INTERVAL_DIV")
  done
  POS_MEAS_COV=$(float_eval "$POS_MEAS_COV / $INTERVAL_DIV")
done

echo "This will perform $MC_TRIALS_COUNT trials"

rm temp_sat3D_airship_Pa_emd.rkx
rm temp_sat3D_airship_R_gyro.rkx
rm $1/tmp_predstats.ssv

else

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


# Command to invoke Octave:
# octave --quiet --norc --no-window-system --eval "add_predstats('$OUTPUTFILENAME','predstats_accum.ssv')" 

fi



