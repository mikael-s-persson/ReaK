#!/bin/bash

echo "Making monte-carlo runs ..."

OUTPUT_DIRECTORY="pp_results/hidim"

BASE_MC_OPTIONS="--mc-runs=100 --mc-results=1 --mc-prog-interval=100 --mc-vertices="
for DIMENSIONS in 3 4 5 6
do
  NUM_VERTICES=1
  for I in $(eval echo {1..$DIMENSIONS});
  do
    ((NUM_VERTICES *= 10 ))
  done
  let "NUM_VERTICES = (NUM_VERTICES > 100000000) ? 100000000 : NUM_VERTICES"
  SPACE_NAME="e$DIMENSIONS"
  MC_OPTIONS=$BASE_MC_OPTIONS$NUM_VERTICES" --mc-space="$SPACE_NAME
  echo "${MC_OPTIONS}"
  
  for STORAGE_I in 1 2
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
      ./test_hidim_planners -o ${OUTPUT_DIRECTORY} $MC_OPTIONS --rrt-star --knn-method=$KNN_METHOD --mg-storage=$STORAGE
      mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_times.txt ${OUTPUT_DIRECTORY}/rrt_star_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_times.txt
      echo "$MC_OPTIONS --rrt-star --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > ${OUTPUT_DIRECTORY}/rrt_star_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_info.txt
#       ./test_hidim_planners -o ${OUTPUT_DIRECTORY} $MC_OPTIONS --rrt-star --rrt-star-with-bnb --knn-method=$KNN_METHOD --mg-storage=$STORAGE
#       mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_times.txt ${OUTPUT_DIRECTORY}/rrt_star_bnb_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_times.txt
#       echo "$MC_OPTIONS --rrt-star --rrt-star-with-bnb --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > ${OUTPUT_DIRECTORY}/$1/rrt_star_bnb_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_info.txt
#       ./test_hidim_planners -o ${OUTPUT_DIRECTORY} $MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --knn-method=$KNN_METHOD --mg-storage=$STORAGE
#       mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_times.txt ${OUTPUT_DIRECTORY}/sbastar_lazy_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_times.txt
#       echo "$MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > ${OUTPUT_DIRECTORY}/$1/sbastar_lazy_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_info.txt
      ./test_hidim_planners -o ${OUTPUT_DIRECTORY} $MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --sba-relaxation=5.0 --knn-method=$KNN_METHOD --mg-storage=$STORAGE
      mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_times.txt ${OUTPUT_DIRECTORY}/sbastar_lazy_any_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_times.txt
      echo "$MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --sba-relaxation=5.0 --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > ${OUTPUT_DIRECTORY}/$1/sbastar_lazy_any_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_info.txt
#       ./test_hidim_planners -o ${OUTPUT_DIRECTORY} $MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --sba-relaxation=5.0 --sba-with-bnb --knn-method=$KNN_METHOD --mg-storage=$STORAGE
#       mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_times.txt ${OUTPUT_DIRECTORY}/sbastar_lazy_any_bnb_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_times.txt
#       echo "$MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --sba-relaxation=5.0 --sba-with-bnb --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > ${OUTPUT_DIRECTORY}/$1/sbastar_lazy_any_bnb_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_info.txt
      ./test_hidim_planners -o ${OUTPUT_DIRECTORY} $MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --sba-relaxation=5.0 --sba-with-voronoi-pull --sba-sa-temperature=2.0 --knn-method=$KNN_METHOD --mg-storage=$STORAGE
      mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_times.txt ${OUTPUT_DIRECTORY}/sbastar_lazy_any_sa_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_times.txt
      echo "$MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --sba-relaxation=5.0 --sba-with-voronoi-pull --sba-sa-temperature=2.0 --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > ${OUTPUT_DIRECTORY}/$1/sbastar_lazy_any_sa_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_info.txt
#       ./test_hidim_planners -o ${OUTPUT_DIRECTORY} $MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --sba-relaxation=5.0 --sba-with-voronoi-pull --sba-sa-temperature=2.0 --sba-with-bnb --knn-method=$KNN_METHOD --mg-storage=$STORAGE
#       mv ${OUTPUT_DIRECTORY}/${SPACE_NAME}_times.txt ${OUTPUT_DIRECTORY}/sbastar_lazy_any_sa_bnb_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_times.txt
#       echo "$MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --sba-relaxation=5.0 --sba-with-voronoi-pull --sba-sa-temperature=2.0 --sba-with-bnb --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > ${OUTPUT_DIRECTORY}/$1/sbastar_lazy_any_sa_bnb_${SPACE_NAME}_${STORAGE_SHORT}${KNN_METHOD}_info.txt
    done
  done
done
