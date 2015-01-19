#!/bin/bash

echo_and_run() { echo "$@" ; "$@" ; }

START_CONFIG="--start-configuration models/CRS_static_start1.rkx"

TARGET_CONFIG="--target-trajectory sim_results/sat3D_airship/$1_traj.rkx"

# SPACE_CONFIG="--space-definition models/space_defs/planning_space_0l.rkx"
# SPACE_CONFIG="--space-definition models/space_defs/planning_space_1c.rkx"
# SPACE_CONFIG="--space-definition models/space_defs/planning_space_2q.rkx"
# SPACE_CONFIG="--space-definition models/space_defs/planning_space_1v.rkx"
# SPACE_CONFIG="--space-definition models/space_defs/planning_space_2a.rkx"
# For temporal SAP space:
SPACE_CONFIG="--space-order 2 --interpolation-method sap --use-temporal-space --use-rate-limited-space --min-travel-dist 0.1 --max-travel-dist 1 --output-space-order 1"

MODEL_CONFIG="--chaser-target-env models/CRS_airship_lab_safer.rkx"

STOP_CRITERIAS="--max-vertices 60 --max-results 50 --prog-interval 1"

BASIC_PLANNER_OPT="--knn-method bf4 --start-delay 10 --max-random-walk 1"

RUNTYPE_CONFIG="--monte-carlo --mc-runs 10"
# RUNTYPE_CONFIG="--single-run"

OUTPUT_DEST_CONFIG="--output-path pp_results/CRS_dynamic --result-file-prefix $1"

PLANNER_CONFIG="--planner-alg rrt_star ${BASIC_PLANNER_OPT}"
echo_and_run ./run_CRS_planner ${OUTPUT_DEST_CONFIG} ${RUNTYPE_CONFIG} ${PLANNER_CONFIG} ${MODEL_CONFIG} ${START_CONFIG} ${TARGET_CONFIG} ${SPACE_CONFIG} ${STOP_CRITERIAS}

PLANNER_CONFIG="--planner-alg sba_star ${BASIC_PLANNER_OPT} --relaxation-factor 5"
echo_and_run ./run_CRS_planner ${OUTPUT_DEST_CONFIG} ${RUNTYPE_CONFIG} ${PLANNER_CONFIG} ${MODEL_CONFIG} ${START_CONFIG} ${TARGET_CONFIG} ${SPACE_CONFIG} ${STOP_CRITERIAS}
mv pp_results/CRS_dynamic/$1/sba_star__any_adj_list_bf4_solutions.txt pp_results/CRS_dynamic/$1/sba_star__any_adj_list_bf4_r05_solutions.txt
mv pp_results/CRS_dynamic/$1/sba_star__any_adj_list_bf4_times.txt pp_results/CRS_dynamic/$1/sba_star__any_adj_list_bf4_r05_times.txt

PLANNER_CONFIG="--planner-alg sba_star ${BASIC_PLANNER_OPT} --relaxation-factor 5 --with-voronoi-pull --sa-temperature 5"
echo_and_run ./run_CRS_planner ${OUTPUT_DEST_CONFIG} ${RUNTYPE_CONFIG} ${PLANNER_CONFIG} ${MODEL_CONFIG} ${START_CONFIG} ${TARGET_CONFIG} ${SPACE_CONFIG} ${STOP_CRITERIAS}
mv pp_results/CRS_dynamic/$1/sba_star__any_sa_adj_list_bf4_solutions.txt pp_results/CRS_dynamic/$1/sba_star__any_sa_adj_list_bf4_r05_t05_solutions.txt
mv pp_results/CRS_dynamic/$1/sba_star__any_sa_adj_list_bf4_times.txt pp_results/CRS_dynamic/$1/sba_star__any_sa_adj_list_bf4_r05_t05_times.txt





