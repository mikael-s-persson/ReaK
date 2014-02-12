#!/bin/bash

echo_and_run() { echo "$@" ; "$@" ; }

START_CONFIG="--start-configuration models/CRS_static_start1.rkx"
TARGET_CONFIG="--target-pose models/CRS_static_target1.rkx"

SPACE_CONFIG="--space-definition models/space_defs/planning_space_1c.rkx"

MODEL_CONFIG="--chaser-target-env models/CRS_static_base_models.pbuf"

STOP_CRITERIAS="--max-vertices 5000 --max-results 50 --prog-interval 10"

RUNTYPE_CONFIG="--monte-carlo --mc-runs 1000"
# RUNTYPE_CONFIG="--single-run"

OUTPUT_DEST_CONFIG="--output-path pp_results/CRS_static --result-file-prefix scene1"

PLANNER_CONFIG="--planner-alg rrt_star --knn-method bf4"
echo_and_run ./run_CRS_planner ${OUTPUT_DEST_CONFIG} ${RUNTYPE_CONFIG} ${PLANNER_CONFIG} ${MODEL_CONFIG} ${START_CONFIG} ${TARGET_CONFIG} ${SPACE_CONFIG} ${STOP_CRITERIAS}

# PLANNER_CONFIG="--planner-alg rrt_star --knn-method bf4 --with-bnb"
# echo_and_run ./run_CRS_planner ${OUTPUT_DEST_CONFIG} ${RUNTYPE_CONFIG} ${PLANNER_CONFIG} ${MODEL_CONFIG} ${START_CONFIG} ${TARGET_CONFIG} ${SPACE_CONFIG} ${STOP_CRITERIAS}

# PLANNER_CONFIG="--planner-alg sba_star --knn-method bf4"
# echo_and_run ./run_CRS_planner ${OUTPUT_DEST_CONFIG} ${RUNTYPE_CONFIG} ${PLANNER_CONFIG} ${MODEL_CONFIG} ${START_CONFIG} ${TARGET_CONFIG} ${SPACE_CONFIG} ${STOP_CRITERIAS}

# PLANNER_CONFIG="--planner-alg sba_star --knn-method bf4 --with-bnb"
# echo_and_run ./run_CRS_planner ${OUTPUT_DEST_CONFIG} ${RUNTYPE_CONFIG} ${PLANNER_CONFIG} ${MODEL_CONFIG} ${START_CONFIG} ${TARGET_CONFIG} ${SPACE_CONFIG} ${STOP_CRITERIAS}

# PLANNER_CONFIG="--planner-alg sba_star --knn-method bf4 --relaxation-factor 10"
# echo_and_run ./run_CRS_planner ${OUTPUT_DEST_CONFIG} ${RUNTYPE_CONFIG} ${PLANNER_CONFIG} ${MODEL_CONFIG} ${START_CONFIG} ${TARGET_CONFIG} ${SPACE_CONFIG} ${STOP_CRITERIAS}

# PLANNER_CONFIG="--planner-alg sba_star --knn-method bf4 --relaxation-factor 10 --with-bnb"
# echo_and_run ./run_CRS_planner ${OUTPUT_DEST_CONFIG} ${RUNTYPE_CONFIG} ${PLANNER_CONFIG} ${MODEL_CONFIG} ${START_CONFIG} ${TARGET_CONFIG} ${SPACE_CONFIG} ${STOP_CRITERIAS}

PLANNER_CONFIG="--planner-alg sba_star --knn-method bf4 --relaxation-factor 10 --with-voronoi-pull --sa-temperature 50"
echo_and_run ./run_CRS_planner ${OUTPUT_DEST_CONFIG} ${RUNTYPE_CONFIG} ${PLANNER_CONFIG} ${MODEL_CONFIG} ${START_CONFIG} ${TARGET_CONFIG} ${SPACE_CONFIG} ${STOP_CRITERIAS}

# PLANNER_CONFIG="--planner-alg sba_star --knn-method bf4 --relaxation-factor 10 --with-voronoi-pull --sa-temperature 50 --with-bnb"
# echo_and_run ./run_CRS_planner ${OUTPUT_DEST_CONFIG} ${RUNTYPE_CONFIG} ${PLANNER_CONFIG} ${MODEL_CONFIG} ${START_CONFIG} ${TARGET_CONFIG} ${SPACE_CONFIG} ${STOP_CRITERIAS}

