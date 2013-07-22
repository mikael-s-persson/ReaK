#!/bin/bash


OUTPUT_DIRECTORY="pp_results/hidim"

if [ $1 -eq 1 ] ;
then
  
  echo "Making single-runs ..."
  
  BASE_SR_OPTIONS="-s --mc-results=50 --mc-prog-interval=100 --mc-vertices="
  SR_KNN_METHOD="bf4"
  SR_STORAGE="adj-list"
  SR_SBA_RELAX="5.0"
  SR_SBA_CUTOFF="0.01"
  for DIMENSIONS in 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
  do
    SR_NUM_VERTICES=100
    for I in $(eval echo {1..$DIMENSIONS});
    do
      ((SR_NUM_VERTICES *= 2 ))
    done
    SR_SA_TEMP=`echo "l( $SR_NUM_VERTICES / 2.0 )" | bc -l`
    let "SR_NUM_VERTICES = (SR_NUM_VERTICES > 1000000) ? 1000000 : SR_NUM_VERTICES"
    SR_SPACE_NAME="e$DIMENSIONS"
    SR_OPTIONS=$BASE_SR_OPTIONS$SR_NUM_VERTICES
    echo "${SR_OPTIONS}"
    LAUNCH_SR_PLANNER_CMD="./test_hidim_planners_"$SR_SPACE_NAME
    
    
    ${LAUNCH_SR_PLANNER_CMD} -o ${OUTPUT_DIRECTORY} $SR_OPTIONS --rrt-star --knn-method=$SR_KNN_METHOD --mg-storage=$SR_STORAGE
    echo "$SR_OPTIONS --rrt-star --knn-method=$SR_KNN_METHOD --mg-storage=$SR_STORAGE" > ${OUTPUT_DIRECTORY}/rrt_star/${SR_SPACE_NAME}_info.txt
    
    ${LAUNCH_SR_PLANNER_CMD} -o ${OUTPUT_DIRECTORY} $SR_OPTIONS --rrt-star --rrt-star-with-bnb --knn-method=$SR_KNN_METHOD --mg-storage=$SR_STORAGE
    echo "$SR_OPTIONS --rrt-star --rrt-star-with-bnb --knn-method=$SR_KNN_METHOD --mg-storage=$SR_STORAGE" > ${OUTPUT_DIRECTORY}/rrt_star_bnb/${SR_SPACE_NAME}_info.txt
    
    ${LAUNCH_SR_PLANNER_CMD} -o ${OUTPUT_DIRECTORY} $SR_OPTIONS --sba-star --sba-density-cutoff=$SR_SBA_CUTOFF --knn-method=$SR_KNN_METHOD --mg-storage=$SR_STORAGE
    echo "$SR_OPTIONS --sba-star --sba-density-cutoff=$SR_SBA_CUTOFF --knn-method=$SR_KNN_METHOD --mg-storage=$SR_STORAGE" > ${OUTPUT_DIRECTORY}/sbastar_lazy/${SR_SPACE_NAME}_info.txt
    
    ${LAUNCH_SR_PLANNER_CMD} -o ${OUTPUT_DIRECTORY} $SR_OPTIONS --sba-star --sba-density-cutoff=$SR_SBA_CUTOFF --sba-relaxation=$SR_SBA_RELAX --knn-method=$SR_KNN_METHOD --mg-storage=$SR_STORAGE
    echo "$SR_OPTIONS --sba-star --sba-density-cutoff=$SR_SBA_CUTOFF --sba-relaxation=$SR_SBA_RELAX --knn-method=$SR_KNN_METHOD --mg-storage=$SR_STORAGE" > ${OUTPUT_DIRECTORY}/sbastar_lazy_any/${SR_SPACE_NAME}_info.txt
    
    ${LAUNCH_SR_PLANNER_CMD} -o ${OUTPUT_DIRECTORY} $SR_OPTIONS --sba-star --sba-density-cutoff=$SR_SBA_CUTOFF --sba-relaxation=$SR_SBA_RELAX --sba-with-bnb --knn-method=$SR_KNN_METHOD --mg-storage=$SR_STORAGE
    echo "$SR_OPTIONS --sba-star --sba-density-cutoff=$SR_SBA_CUTOFF --sba-relaxation=$SR_SBA_RELAX --sba-with-bnb --knn-method=$SR_KNN_METHOD --mg-storage=$SR_STORAGE" > ${OUTPUT_DIRECTORY}/sbastar_lazy_any_bnb/${SR_SPACE_NAME}_info.txt
    
    ${LAUNCH_SR_PLANNER_CMD} -o ${OUTPUT_DIRECTORY} $SR_OPTIONS --sba-star --sba-density-cutoff=$SR_SBA_CUTOFF --sba-relaxation=$SR_SBA_RELAX --sba-with-voronoi-pull --sba-sa-temperature=$SR_SA_TEMP --knn-method=$SR_KNN_METHOD --mg-storage=$SR_STORAGE
    echo "$SR_OPTIONS --sba-star --sba-density-cutoff=$SR_SBA_CUTOFF --sba-relaxation=$SR_SBA_RELAX --sba-with-voronoi-pull --sba-sa-temperature=$SR_SA_TEMP --knn-method=$SR_KNN_METHOD --mg-storage=$SR_STORAGE" > ${OUTPUT_DIRECTORY}/sbastar_lazy_any_sa/${SR_SPACE_NAME}_info.txt
    
    ${LAUNCH_SR_PLANNER_CMD} -o ${OUTPUT_DIRECTORY} $SR_OPTIONS --sba-star --sba-density-cutoff=$SR_SBA_CUTOFF --sba-relaxation=$SR_SBA_RELAX --sba-with-voronoi-pull --sba-sa-temperature=$SR_SA_TEMP --sba-with-bnb --knn-method=$SR_KNN_METHOD --mg-storage=$SR_STORAGE
    echo "$SR_OPTIONS --sba-star --sba-density-cutoff=$SR_SBA_CUTOFF --sba-relaxation=$SR_SBA_RELAX --sba-with-voronoi-pull --sba-sa-temperature=$SR_SA_TEMP --sba-with-bnb --knn-method=$SR_KNN_METHOD --mg-storage=$SR_STORAGE" > ${OUTPUT_DIRECTORY}/sbastar_lazy_any_sa_bnb/${SR_SPACE_NAME}_info.txt
    
  done
  
fi


if [ $1 -eq 2 ] ;
then
  
  echo "Making monte-carlo runs ..."
  
  BASE_MC_OPTIONS="-m --mc-runs=100 --mc-results=1 --mc-prog-interval=10 --mc-vertices="
  MC_SBA_RELAX="5.0"
  MC_SBA_CUTOFF="0.01"
  for DIMENSIONS in 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
  do
    NUM_VERTICES=100
    for I in $(eval echo {1..$DIMENSIONS});
    do
      ((NUM_VERTICES *= 2 ))
    done
    MC_SA_TEMP=`echo "l( $NUM_VERTICES / 2.0 )" | bc -l`
    echo "${MC_SA_TEMP}"
    let "NUM_VERTICES = (NUM_VERTICES > 1000000) ? 1000000 : NUM_VERTICES"
    SPACE_NAME="e$DIMENSIONS"
  #   MC_OPTIONS=$BASE_MC_OPTIONS$NUM_VERTICES" --mc-space="$SPACE_NAME
    MC_OPTIONS=$BASE_MC_OPTIONS$NUM_VERTICES
    echo "${MC_OPTIONS}"
    LAUNCH_PLANNER_CMD="./test_hidim_planners_"$SPACE_NAME
    
    for STORAGE_I in 1
    do
      STORAGE="adj-list"
      STORAGE_SHORT="a"
      if [ $STORAGE_I -eq 2 ] ;
      then 
        STORAGE="dvp-adj-list"
        STORAGE_SHORT="d"
      fi
      for KNN_METHOD in bf4
      do

        ${LAUNCH_PLANNER_CMD} -o ${OUTPUT_DIRECTORY} $MC_OPTIONS --rrt-star --knn-method=$KNN_METHOD --mg-storage=$STORAGE
        mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_times.txt ${OUTPUT_DIRECTORY}/rrt_star_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_times.txt
        mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_solutions.txt ${OUTPUT_DIRECTORY}/rrt_star_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
        echo "$MC_OPTIONS --rrt-star --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > ${OUTPUT_DIRECTORY}/rrt_star_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_info.txt

  #       ${LAUNCH_PLANNER_CMD} -o ${OUTPUT_DIRECTORY} $MC_OPTIONS --rrt-star --rrt-star-with-bnb --knn-method=$KNN_METHOD --mg-storage=$STORAGE
  #       mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_times.txt ${OUTPUT_DIRECTORY}/rrt_star_bnb_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_times.txt
  #       mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_solutions.txt ${OUTPUT_DIRECTORY}/rrt_star_bnb_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
  #       echo "$MC_OPTIONS --rrt-star --rrt-star-with-bnb --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > ${OUTPUT_DIRECTORY}/rrt_star_bnb_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_info.txt

  #       ${LAUNCH_PLANNER_CMD} -o ${OUTPUT_DIRECTORY} $MC_OPTIONS --sba-star --sba-density-cutoff=$MC_SBA_CUTOFF --knn-method=$KNN_METHOD --mg-storage=$STORAGE
  #       mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_times.txt ${OUTPUT_DIRECTORY}/sbastar_lazy_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_times.txt
  #       mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_solutions.txt ${OUTPUT_DIRECTORY}/sbastar_lazy_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
  #       echo "$MC_OPTIONS --sba-star --sba-density-cutoff=$MC_SBA_CUTOFF --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > ${OUTPUT_DIRECTORY}/sbastar_lazy_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_info.txt

#        ${LAUNCH_PLANNER_CMD} -o ${OUTPUT_DIRECTORY} $MC_OPTIONS --sba-star --sba-density-cutoff=$MC_SBA_CUTOFF --sba-relaxation=$MC_SBA_RELAX --knn-method=$KNN_METHOD --mg-storage=$STORAGE
#        mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_times.txt ${OUTPUT_DIRECTORY}/sbastar_lazy_any_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_times.txt
#        mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_solutions.txt ${OUTPUT_DIRECTORY}/sbastar_lazy_any_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
#        echo "$MC_OPTIONS --sba-star --sba-density-cutoff=$MC_SBA_CUTOFF --sba-relaxation=$MC_SBA_RELAX --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > ${OUTPUT_DIRECTORY}/sbastar_lazy_any_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_info.txt

  #       ${LAUNCH_PLANNER_CMD} -o ${OUTPUT_DIRECTORY} $MC_OPTIONS --sba-star --sba-density-cutoff=$MC_SBA_CUTOFF --sba-relaxation=$MC_SBA_RELAX --sba-with-bnb --knn-method=$KNN_METHOD --mg-storage=$STORAGE
  #       mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_times.txt ${OUTPUT_DIRECTORY}/sbastar_lazy_any_bnb_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_times.txt
  #       mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_solutions.txt ${OUTPUT_DIRECTORY}/sbastar_lazy_any_bnb_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
  #       echo "$MC_OPTIONS --sba-star --sba-density-cutoff=$MC_SBA_CUTOFF --sba-relaxation=$MC_SBA_RELAX --sba-with-bnb --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > ${OUTPUT_DIRECTORY}/sbastar_lazy_any_bnb_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_info.txt

#        ${LAUNCH_PLANNER_CMD} -o ${OUTPUT_DIRECTORY} $MC_OPTIONS --sba-star --sba-density-cutoff=$MC_SBA_CUTOFF --sba-relaxation=$MC_SBA_RELAX --sba-with-voronoi-pull --sba-sa-temperature=$MC_SA_TEMP --knn-method=$KNN_METHOD --mg-storage=$STORAGE
#        mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_times.txt ${OUTPUT_DIRECTORY}/sbastar_lazy_any_sa_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_times.txt
#        mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_solutions.txt ${OUTPUT_DIRECTORY}/sbastar_lazy_any_sa_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
#        echo "$MC_OPTIONS --sba-star --sba-density-cutoff=$MC_SBA_CUTOFF --sba-relaxation=$MC_SBA_RELAX --sba-with-voronoi-pull --sba-sa-temperature=$MC_SA_TEMP --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > ${OUTPUT_DIRECTORY}/sbastar_lazy_any_sa_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_info.txt

  #       ${LAUNCH_PLANNER_CMD} -o ${OUTPUT_DIRECTORY} $MC_OPTIONS --sba-star --sba-density-cutoff=$MC_SBA_CUTOFF --sba-relaxation=$MC_SBA_RELAX --sba-with-voronoi-pull --sba-sa-temperature=$MC_SA_TEMP --sba-with-bnb --knn-method=$KNN_METHOD --mg-storage=$STORAGE
  #       mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_times.txt ${OUTPUT_DIRECTORY}/sbastar_lazy_any_sa_bnb_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_times.txt
  #       mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_solutions.txt ${OUTPUT_DIRECTORY}/sbastar_lazy_any_sa_bnb_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
  #       echo "$MC_OPTIONS --sba-star --sba-density-cutoff=$MC_SBA_CUTOFF --sba-relaxation=$MC_SBA_RELAX --sba-with-voronoi-pull --sba-sa-temperature=$MC_SA_TEMP --sba-with-bnb --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > ${OUTPUT_DIRECTORY}/sbastar_lazy_any_sa_bnb_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_info.txt

      done
    done
  done
fi