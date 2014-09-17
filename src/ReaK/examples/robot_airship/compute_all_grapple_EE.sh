#!/bin/bash

for f in $1meas_*.ssv
do
  exp_num=${f//meas_}
  exp_num=${exp_num//.ssv}
  exp_num=${exp_num//$1}
  echo "Processing target record ${exp_num}..."
  ./pose_to_grapple --input $f --output ${f//meas_/grapple_}
  ./pose_to_grapple --input ${f//meas_/pred_} --output ${f//meas_/grapple_pred_}
  ./pose_to_grapple --input ${f//meas_/est_} --output ${f//meas_/grapple_est_}
done

for f in $1jtctrl_*.ssv
do
  exp_num=${f//jtctrl_}
  exp_num=${exp_num//.ssv}
  exp_num=${exp_num//$1}
  echo "Processing jtctrl record ${exp_num}..."
  ./CRS_jtctrl_to_EE --input $f --output ${f//jtctrl_/EE_}
done

