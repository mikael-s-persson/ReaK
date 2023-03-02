#!/ bin / bash

echo_and_run() {
  echo "$@";
  "$@";
}

START_CONFIG =
    "--start-configuration models/CRS_static_start1.rkx" TARGET_CONFIG =
        "--target-pose models/CRS_static_target2.rkx"

#SPACE_CONFIG = "--space-definition models/space_defs/planning_space_0l.rkx"
    SPACE_CONFIG = "--space-definition models/space_defs/planning_space_1c.rkx"
#SPACE_CONFIG = "--space-definition models/space_defs/planning_space_2q.rkx"
#SPACE_CONFIG = "--space-definition models/space_defs/planning_space_1v.rkx"
#SPACE_CONFIG = "--space-definition models/space_defs/planning_space_2a.rkx"

    MODEL_CONFIG = "--chaser-target-env models/CRS_static_base_models.pbuf"

    STOP_CRITERIAS = "--max-vertices 6000 --max-results 50 --prog-interval 1"

    RUNTYPE_CONFIG = "--monte-carlo --mc-runs 1000"
#RUNTYPE_CONFIG = "--single-run"

    OUTPUT_DEST_CONFIG =
        "--output-path pp_results/CRS_static --result-file-prefix scene2_1c"

    PLANNER_CONFIG =
        "--planner-alg rrt_star --knn-method bf4" echo_and_run./
        run_CRS_planner ${OUTPUT_DEST_CONFIG} ${RUNTYPE_CONFIG} ${
            PLANNER_CONFIG} ${MODEL_CONFIG} ${START_CONFIG} ${TARGET_CONFIG} ${
            SPACE_CONFIG} ${STOP_CRITERIAS}

        PLANNER_CONFIG =
            "--planner-alg sba_star --knn-method bf4 --relaxation-factor "
            "5" echo_and_run./
            run_CRS_planner ${OUTPUT_DEST_CONFIG} ${RUNTYPE_CONFIG} ${
                PLANNER_CONFIG} ${MODEL_CONFIG} ${START_CONFIG} ${
                TARGET_CONFIG} ${SPACE_CONFIG} ${STOP_CRITERIAS} mv pp_results /
            CRS_static / scene2_1c /
            sba_star__any_adj_list_bf4_solutions.txt pp_results / CRS_static /
            scene2_1c /
            sba_star__any_adj_list_bf4_r05_solutions.txt mv pp_results /
            CRS_static / scene2_1c /
            sba_star__any_adj_list_bf4_times.txt pp_results / CRS_static /
            scene2_1c /
            sba_star__any_adj_list_bf4_r05_times.txt

                PLANNER_CONFIG =
                "--planner-alg sba_star --knn-method bf4 --relaxation-factor "
                "10" echo_and_run./
                run_CRS_planner ${OUTPUT_DEST_CONFIG} ${RUNTYPE_CONFIG} ${
                    PLANNER_CONFIG} ${MODEL_CONFIG} ${START_CONFIG} ${
                    TARGET_CONFIG} ${SPACE_CONFIG} ${
                    STOP_CRITERIAS} mv pp_results /
                CRS_static / scene2_1c /
                sba_star__any_adj_list_bf4_solutions.txt pp_results /
                CRS_static / scene2_1c /
                sba_star__any_adj_list_bf4_r10_solutions.txt mv pp_results /
                CRS_static / scene2_1c /
                sba_star__any_adj_list_bf4_times.txt pp_results / CRS_static /
                scene2_1c / sba_star__any_adj_list_bf4_r10_times.txt
