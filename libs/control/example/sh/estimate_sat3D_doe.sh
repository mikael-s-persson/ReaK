#!/ bin / bash

#All the constant options here:

OPTIONS="--init-motion models/satellite3D_init.rkx"
OPTIONS="${OPTIONS} --inertia models/sat3D_airship_inertia.rkx"
OPTIONS="${OPTIONS} --Q-matrix models/sat3D_airship_Q.rkx"
OPTIONS="${OPTIONS} --IMU-config models/sat3D_airship_IMU_config.rkx"
OPTIONS="${OPTIONS} --time-step 0.01"
#OPTIONS = "${OPTIONS} --imkf"
#OPTIONS = "${OPTIONS} --imkf-em"
OPTIONS="${OPTIONS} --imkf-emd"
#OPTIONS = "${OPTIONS} --imkf-emdJ"

#OPTIONS = "${OPTIONS} --R-matrix temp_sat3D_airship_R.rkx"

OPTIONS="${OPTIONS} --R-matrix temp_sat3D_airship_R_gyro.rkx"
OPTIONS="${OPTIONS} --gyro"

OPTIONS="${OPTIONS} --prediction-runs --prediction-interval 5.0"

OPTIONS="${OPTIONS} --pred-assumption 0"
#OPTIONS = "${OPTIONS} --pred-assumption 1"

#OPTIONS = "${OPTIONS} --tsosakf --Pa-matrix temp_sat3D_airship_Pa_em.rkx"
OPTIONS="${OPTIONS} --tsosakf --Pa-matrix temp_sat3D_airship_Pa_emd.rkx"
#OPTIONS = "${OPTIONS} --tsosakf --Pa-matrix temp_sat3D_airship_Pa_emdJ.rkx"

OPTIONS="${OPTIONS} --monte-carlo --min-skips 3 --max-skips 5"

if [[ "0" == "0" ]]
then

INTERVAL_DIV=4.0

DM_COV_LO=0.00001
DM_COV_HI=0.001

ECC_COV_LO=0.0000009
ECC_COV_HI=0.000001

LD_COV_LO=10
LD_COV_HI=1000

RD_COV_LO=10
RD_COV_HI=1000

POS_MEAS_COV_LO=0.000005
POS_MEAS_COV_HI=0.0005
ANG_MEAS_COV_LO=0.000009
ANG_MEAS_COV_HI=0.0009
GYRO_MEAS_COV_LO=0.000005
GYRO_MEAS_COV_HI=0.0005

MEAS_COV_FACT_LO=0.9
MEAS_COV_FACT_HI=1

FINAL_OUTPUT_FILE="$1/pred/all_predstats.ssv"

# #Start from the beginning:
#DM_COV = $DM_COV_HI
#ECC_COV = $ECC_COV_HI
#LD_COV = $LD_COV_HI
#RD_COV = $RD_COV_HI
#
#POS_MEAS_COV = $POS_MEAS_COV_HI
#ANG_MEAS_COV = $ANG_MEAS_COV_HI
#GYRO_MEAS_COV = $GYRO_MEAS_COV_HI
#MEAS_COV_FACT = $MEAS_COV_FACT_HI
#
#echo                                                                                                                                                               \
    "% pos_meas_cov ang_meas_cov gyro_meas_cov dm_cov ecc_cov ld_cov rd_cov P_th success_rate pred_start_time ep_m ea_m ev_m ew_m pdf_est lr_est pdf_meas lr_meas"> \
    $FINAL_OUTPUT_FILE

#Start from this point:
POS_MEAS_COV=0.0005
ANG_MEAS_COV=0.0009
GYRO_MEAS_COV=0.0005
MEAS_COV_FACT=1

DM_COV=.00001562500000000000
ECC_COV=0.000001
LD_COV=250.00000000000000000000
RD_COV=62.50000000000000000000


mkdir $1/pred

(( MC_TRIALS_COUNT = 0 ))


until (( $(echo "$MEAS_COV_FACT > $MEAS_COV_FACT_LO" | bc -q 2>/dev/null) == "0" ))
do
  echo "meas cov factor = $MEAS_COV_FACT"

#until(($(echo "$POS_MEAS_COV > $POS_MEAS_COV_LO" | \
          bc - q 2 > / dev / null) == "0"))
#do
#echo "pos-meas cov = $POS_MEAS_COV"
#until(($(echo "$ANG_MEAS_COV > $ANG_MEAS_COV_LO" | \
          bc - q 2 > / dev / null) == "0"))
#do
#echo "ang-meas cov = $ANG_MEAS_COV"
#until(($(echo "$GYRO_MEAS_COV > $GYRO_MEAS_COV_LO" | \
          bc - q 2 > / dev / null) == "0"))
#do
#echo "gyro-meas cov = $GYRO_MEAS_COV"

      until (( $(echo "$DM_COV > $DM_COV_LO" | bc -q 2>/dev/null) == "0" ))
      do
        echo "dm cov = $DM_COV"
        until (( $(echo "$ECC_COV > $ECC_COV_LO" | bc -q 2>/dev/null) == "0" ))
        do
          echo "ecc. cov = $ECC_COV"
          until (( $(echo "$LD_COV > $LD_COV_LO" | bc -q 2>/dev/null) == "0" ))
          do
            echo "lin-drag cov = $LD_COV"
#RD_COV = $LD_COV
            until (( $(echo "$RD_COV > $RD_COV_LO" | bc -q 2>/dev/null) == "0" ))
            do
              echo "rot-drag cov = $RD_COV"
              
              if [[ "0" == "0" ]]
              then
              ACCUM_FILE="$1/pred/predstats_${POS_MEAS_COV}_${ANG_MEAS_COV}_${GYRO_MEAS_COV}_${DM_COV}_${ECC_COV}_${LD_COV}_${RD_COV}.ssv"
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
              
              for f in $1/*.ssv
              do
                echo "Processing file: $f"
                IN_OPTIONS="--input $f --input-format ssv"
                OUT_OPTIONS="--output $1/pred/tmp_predstats.ssv --output-format ssv"
                ./estimate_satellite3D ${OPTIONS} ${IN_OPTIONS} ${OUT_OPTIONS}
                for g in $1/pred/tmp_predstats*.ssv
                do
                  octave --quiet --norc --no-window-system --eval "add_predstats('$g','$ACCUM_FILE')"
                done
              done
              
              # Put all accumulated data into a single file:
              while read p; do
                echo "$POS_MEAS_COV $ANG_MEAS_COV $GYRO_MEAS_COV $DM_COV $ECC_COV $LD_COV $RD_COV $p" >> $FINAL_OUTPUT_FILE
              done <$ACCUM_FILE
              
              else
              
              (( MC_TRIALS_COUNT++ ))
              
              fi
              
              RD_COV=$(echo "scale=20; $RD_COV / $INTERVAL_DIV" | bc -q 2>/dev/null)
            done
            RD_COV=$RD_COV_HI
            LD_COV=$(echo "scale=20; $LD_COV / $INTERVAL_DIV" | bc -q 2>/dev/null)
          done
          LD_COV=$LD_COV_HI
          ECC_COV=$(echo "scale=20; $ECC_COV / $INTERVAL_DIV" | bc -q 2>/dev/null)
        done
        ECC_COV=$ECC_COV_HI
        DM_COV=$(echo "scale=20; $DM_COV / $INTERVAL_DIV" | bc -q 2>/dev/null)
      done
      DM_COV=$DM_COV_HI
#       GYRO_MEAS_COV=$(echo "scale=20; $GYRO_MEAS_COV / $INTERVAL_DIV" | bc -q 2>/dev/null)
#     done
#     GYRO_MEAS_COV=$GYRO_MEAS_COV_HI
#     ANG_MEAS_COV=$(echo "scale=20; $ANG_MEAS_COV / $INTERVAL_DIV" | bc -q 2>/dev/null)
#   done
#   ANG_MEAS_COV=$ANG_MEAS_COV_HI
#   POS_MEAS_COV=$(echo "scale=20; $POS_MEAS_COV / $INTERVAL_DIV" | bc -q 2>/dev/null)
# done
# POS_MEAS_COV=$POS_MEAS_COV_HI

  MEAS_COV_FACT=$(echo "scale=20; $MEAS_COV_FACT / $INTERVAL_DIV" | bc -q 2>/dev/null)
  POS_MEAS_COV=$(echo "scale=20; $MEAS_COV_FACT * 0.00005" | bc -q 2>/dev/null)
  ANG_MEAS_COV=$(echo "scale=20; $MEAS_COV_FACT * 0.00009" | bc -q 2>/dev/null)
  GYRO_MEAS_COV=$(echo "scale=20; $MEAS_COV_FACT * 0.00005" | bc -q 2>/dev/null)
done

echo "This will perform $MC_TRIALS_COUNT trials"

rm temp_sat3D_airship_Pa_emd.rkx
rm temp_sat3D_airship_R_gyro.rkx
rm $1/pred/tmp_predstats*.ssv

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



