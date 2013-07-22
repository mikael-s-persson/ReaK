#!/bin/bash

if [ $2 -eq 1 ] ;
then
  echo "Making single runs ..."
  SR_OPTIONS="-s --max-vertices=3000 --max-results=50 --prog-interval=10"
  SR_STRUCTS="--knn-method=bf4 --mg-storage=adj-list"

  ./test_planners -i $1.bmp -o pp_results/$1 $SR_OPTIONS --rrt $SR_STRUCTS
  echo "$SR_OPTIONS --rrt $SR_STRUCTS" > pp_results/$1/rrt/info.txt

  ./test_planners -i $1.bmp -o pp_results/$1 $SR_OPTIONS --bi-rrt $SR_STRUCTS
  echo "$SR_OPTIONS --bi-rrt $SR_STRUCTS" > pp_results/$1/birrt/info.txt

  ./test_planners -i $1.bmp -o pp_results/$1 $SR_OPTIONS --rrt-star $SR_STRUCTS
  echo "$SR_OPTIONS --rrt-star $SR_STRUCTS" > pp_results/$1/rrt_star/info.txt

  ./test_planners -i $1.bmp -o pp_results/$1 $SR_OPTIONS --rrt-star --rrt-star-with-bnb $SR_STRUCTS
  echo "$SR_OPTIONS --rrt-star --rrt-star-with-bnb $SR_STRUCTS" > pp_results/$1/rrt_star_bnb/info.txt

  ./test_planners -i $1.bmp -o pp_results/$1 $SR_OPTIONS --prm $SR_STRUCTS
  echo "$SR_OPTIONS --prm $SR_STRUCTS" > pp_results/$1/prm/info.txt

#   ./test_planners -i $1.bmp -o pp_results/$1 $SR_OPTIONS --fadprm --fadprm-relaxation=10.0 $SR_STRUCTS
#   echo "$SR_OPTIONS --fadprm --fadprm-relaxation=10.0 $SR_STRUCTS" > pp_results/$1/fadprm/info.txt

  ./test_planners -i $1.bmp -o pp_results/$1 $SR_OPTIONS --sba-star --sba-density-cutoff=0.01 $SR_STRUCTS
  echo "$SR_OPTIONS --sba-star --sba-density-cutoff=0.01 $SR_STRUCTS" > pp_results/$1/sbastar_lazy/info.txt

  ./test_planners -i $1.bmp -o pp_results/$1 $SR_OPTIONS --sba-star --sba-density-cutoff=0.01 --sba-relaxation=5.0 $SR_STRUCTS
  echo "$SR_OPTIONS --sba-star --sba-density-cutoff=0.01 --sba-relaxation=5.0 $SR_STRUCTS" > pp_results/$1/sbastar_lazy_any/info.txt

  ./test_planners -i $1.bmp -o pp_results/$1 $SR_OPTIONS --sba-star --sba-density-cutoff=0.01 --sba-relaxation=5.0 --sba-with-bnb $SR_STRUCTS
  echo "$SR_OPTIONS --sba-star --sba-density-cutoff=0.01 --sba-relaxation=5.0 --sba-with-bnb $SR_STRUCTS" > pp_results/$1/sbastar_lazy_any_bnb/info.txt

  ./test_planners -i $1.bmp -o pp_results/$1 $SR_OPTIONS --sba-star --sba-density-cutoff=0.01 --sba-relaxation=5.0 --sba-with-voronoi-pull $SR_STRUCTS
  mv pp_results/$1/sbastar_lazy_any_sa pp_results/$1/sbastar_rrtstar
  echo "$SR_OPTIONS --sba-star --sba-density-cutoff=0.01 --sba-relaxation=5.0 --sba-with-voronoi-pull $SR_STRUCTS" > pp_results/$1/sbastar_rrtstar/info.txt

  ./test_planners -i $1.bmp -o pp_results/$1 $SR_OPTIONS --sba-star --sba-density-cutoff=0.01 --sba-relaxation=5.0 --sba-with-voronoi-pull --sba-sa-temperature=2.0 $SR_STRUCTS
  echo "$SR_OPTIONS --sba-star --sba-density-cutoff=0.01 --sba-relaxation=5.0 --sba-with-voronoi-pull --sba-sa-temperature=2.0 $SR_STRUCTS" > pp_results/$1/sbastar_lazy_any_sa/info.txt

  ./test_planners -i $1.bmp -o pp_results/$1 $SR_OPTIONS --sba-star --sba-density-cutoff=0.01 --sba-relaxation=5.0 --sba-with-voronoi-pull --sba-sa-temperature=2.0 --sba-with-bnb $SR_STRUCTS
  echo "$SR_OPTIONS --sba-star --sba-density-cutoff=0.01 --sba-relaxation=5.0 --sba-with-voronoi-pull --sba-sa-temperature=2.0 --sba-with-bnb $SR_STRUCTS" > pp_results/$1/sbastar_lazy_any_sa_bnb/info.txt

fi

if [ $2 -eq 2 ] ;
then
  echo "Making monte-carlo runs ..."
  MC_OPTIONS="-m --mc-runs=1000 --mc-vertices=3000 --mc-results=50 --mc-prog-interval=10"
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

      ./test_planners -i $1.bmp -o pp_results/$1 $MC_OPTIONS --rrt --knn-method=$KNN_METHOD --mg-storage=$STORAGE
      mv pp_results/$1/$1_times.txt pp_results/$1/rrt_${STORAGE_SHORT}${KNN_METHOD}_times.txt
      mv pp_results/$1/$1_solutions.txt pp_results/$1/rrt_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
      echo "$MC_OPTIONS --rrt --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > pp_results/$1/rrt_${STORAGE_SHORT}${KNN_METHOD}_info.txt

      ./test_planners -i $1.bmp -o pp_results/$1 $MC_OPTIONS --bi-rrt --knn-method=$KNN_METHOD --mg-storage=$STORAGE
      mv pp_results/$1/$1_times.txt pp_results/$1/birrt_${STORAGE_SHORT}${KNN_METHOD}_times.txt
      mv pp_results/$1/$1_solutions.txt pp_results/$1/birrt_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
      echo "$MC_OPTIONS --bi-rrt --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > pp_results/$1/birrt_${STORAGE_SHORT}${KNN_METHOD}_info.txt

      ./test_planners -i $1.bmp -o pp_results/$1 $MC_OPTIONS --rrt-star --knn-method=$KNN_METHOD --mg-storage=$STORAGE
      mv pp_results/$1/$1_times.txt pp_results/$1/rrt_star_${STORAGE_SHORT}${KNN_METHOD}_times.txt
      mv pp_results/$1/$1_solutions.txt pp_results/$1/rrt_star_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
      echo "$MC_OPTIONS --rrt-star --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > pp_results/$1/rrt_star_${STORAGE_SHORT}${KNN_METHOD}_info.txt

#       ./test_planners -i $1.bmp -o pp_results/$1 $MC_OPTIONS --rrt-star --rrt-star-with-bnb --knn-method=$KNN_METHOD --mg-storage=$STORAGE
#       mv pp_results/$1/$1_times.txt pp_results/$1/rrt_star_bnb_${STORAGE_SHORT}${KNN_METHOD}_times.txt
#       mv pp_results/$1/$1_solutions.txt pp_results/$1/rrt_star_bnb_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
#       echo "$MC_OPTIONS --rrt-star --rrt-star-with-bnb --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > pp_results/$1/rrt_star_bnb_${STORAGE_SHORT}${KNN_METHOD}_info.txt

      ./test_planners -i $1.bmp -o pp_results/$1 $MC_OPTIONS --prm --knn-method=$KNN_METHOD --mg-storage=$STORAGE
      mv pp_results/$1/$1_times.txt pp_results/$1/prm_${STORAGE_SHORT}${KNN_METHOD}_times.txt
      mv pp_results/$1/$1_solutions.txt pp_results/$1/prm_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
      echo "$MC_OPTIONS --prm --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > pp_results/$1/prm_${STORAGE_SHORT}${KNN_METHOD}_info.txt

#       ./test_planners -i $1.bmp -o pp_results/$1 $MC_OPTIONS --fadprm --fadprm-relaxation=10.0 --knn-method=$KNN_METHOD --mg-storage=$STORAGE
#       mv pp_results/$1/$1_times.txt pp_results/$1/fadprm_${STORAGE_SHORT}${KNN_METHOD}_times.txt
#       mv pp_results/$1/$1_solutions.txt pp_results/$1/fadprm_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
#       echo "$MC_OPTIONS --fadprm --fadprm-relaxation=10.0 --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > pp_results/$1/fadprm_${STORAGE_SHORT}${KNN_METHOD}_info.txt

#       ./test_planners -i $1.bmp -o pp_results/$1 $MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --knn-method=$KNN_METHOD --mg-storage=$STORAGE
#       mv pp_results/$1/$1_times.txt pp_results/$1/sbastar_lazy_${STORAGE_SHORT}${KNN_METHOD}_times.txt
#       mv pp_results/$1/$1_solutions.txt pp_results/$1/sbastar_lazy_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
#       echo "$MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > pp_results/$1/sbastar_lazy_${STORAGE_SHORT}${KNN_METHOD}_info.txt

      ./test_planners -i $1.bmp -o pp_results/$1 $MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --sba-relaxation=5.0 --knn-method=$KNN_METHOD --mg-storage=$STORAGE
      mv pp_results/$1/$1_times.txt pp_results/$1/sbastar_lazy_any_${STORAGE_SHORT}${KNN_METHOD}_times.txt
      mv pp_results/$1/$1_solutions.txt pp_results/$1/sbastar_lazy_any_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
      echo "$MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --sba-relaxation=5.0 --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > pp_results/$1/sbastar_lazy_any_${STORAGE_SHORT}${KNN_METHOD}_info.txt

#       ./test_planners -i $1.bmp -o pp_results/$1 $MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --sba-relaxation=5.0 --sba-with-bnb --knn-method=$KNN_METHOD --mg-storage=$STORAGE
#       mv pp_results/$1/$1_times.txt pp_results/$1/sbastar_lazy_any_bnb_${STORAGE_SHORT}${KNN_METHOD}_times.txt
#       mv pp_results/$1/$1_solutions.txt pp_results/$1/sbastar_lazy_any_bnb_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
#       echo "$MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --sba-relaxation=5.0 --sba-with-bnb --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > pp_results/$1/sbastar_lazy_any_bnb_${STORAGE_SHORT}${KNN_METHOD}_info.txt

      ./test_planners -i $1.bmp -o pp_results/$1 $MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --sba-relaxation=5.0 --sba-with-voronoi-pull --knn-method=$KNN_METHOD --mg-storage=$STORAGE
      mv pp_results/$1/$1_times.txt pp_results/$1/sbastar_rrtstar_${STORAGE_SHORT}${KNN_METHOD}_times.txt
      mv pp_results/$1/$1_solutions.txt pp_results/$1/sbastar_rrtstar_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
      echo "$MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --sba-relaxation=5.0 --sba-with-voronoi-pull --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > pp_results/$1/sbastar_rrtstar_${STORAGE_SHORT}${KNN_METHOD}_info.txt

      ./test_planners -i $1.bmp -o pp_results/$1 $MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --sba-relaxation=5.0 --sba-with-voronoi-pull --sba-sa-temperature=2.0 --knn-method=$KNN_METHOD --mg-storage=$STORAGE
      mv pp_results/$1/$1_times.txt pp_results/$1/sbastar_lazy_any_sa_${STORAGE_SHORT}${KNN_METHOD}_times.txt
      mv pp_results/$1/$1_solutions.txt pp_results/$1/sbastar_lazy_any_sa_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
      echo "$MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --sba-relaxation=5.0 --sba-with-voronoi-pull --sba-sa-temperature=2.0 --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > pp_results/$1/sbastar_lazy_any_sa_${STORAGE_SHORT}${KNN_METHOD}_info.txt

#       ./test_planners -i $1.bmp -o pp_results/$1 $MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --sba-relaxation=5.0 --sba-with-voronoi-pull --sba-sa-temperature=2.0 --sba-with-bnb --knn-method=$KNN_METHOD --mg-storage=$STORAGE
#       mv pp_results/$1/$1_times.txt pp_results/$1/sbastar_lazy_any_sa_bnb_${STORAGE_SHORT}${KNN_METHOD}_times.txt
#       mv pp_results/$1/$1_solutions.txt pp_results/$1/sbastar_lazy_any_sa_bnb_${STORAGE_SHORT}${KNN_METHOD}_solutions.txt
#       echo "$MC_OPTIONS --sba-star --sba-potential-cutoff=0.02 --sba-density-cutoff=0.0 --sba-relaxation=5.0 --sba-with-voronoi-pull --sba-sa-temperature=2.0 --sba-with-bnb --knn-method=$KNN_METHOD --mg-storage=$STORAGE" > pp_results/$1/sbastar_lazy_any_sa_bnb_${STORAGE_SHORT}${KNN_METHOD}_info.txt

    done
  done
fi
