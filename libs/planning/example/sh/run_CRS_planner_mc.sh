#!/ bin / bash

echo_and_run() {
  echo "$@";
  "$@";
}

START_CONFIG =
    "--start-configuration models/CRS_static_start1.rkx" TARGET_CONFIG =
        "--target-pose models/CRS_static_target2.rkx"

    SPACE_CONFIG = "--space-definition models/space_defs/planning_space_0l.rkx"
#SPACE_CONFIG = "--space-definition models/space_defs/planning_space_1c.rkx"
#SPACE_CONFIG = "--space-definition models/space_defs/planning_space_2q.rkx"
#SPACE_CONFIG = "--space-definition models/space_defs/planning_space_1v.rkx"
#SPACE_CONFIG = "--space-definition models/space_defs/planning_space_2a.rkx"

    MODEL_CONFIG = "--chaser-target-env models/CRS_static_base_models.pbuf"

    STOP_CRITERIAS = "--max-vertices 3000 --max-results 50 --prog-interval 1"

    RUNTYPE_CONFIG = "--monte-carlo --mc-runs 100"
#RUNTYPE_CONFIG = "--single-run"

    OUTPUT_DEST_CONFIG =
        "--output-path pp_results/CRS_static --result-file-prefix scene2"

    PLANNER_CONFIG =
        "--planner-alg rrt_star --knn-method bf4" echo_and_run./
        run_CRS_planner ${OUTPUT_DEST_CONFIG} ${RUNTYPE_CONFIG} ${
            PLANNER_CONFIG} ${MODEL_CONFIG} ${START_CONFIG} ${TARGET_CONFIG} ${
            SPACE_CONFIG} ${STOP_CRITERIAS}

#PLANNER_CONFIG = "--planner-alg rrt_star --knn-method bf4 --with-bnb"
#echo_and_run./ run_CRS_planner ${OUTPUT_DEST_CONFIG } ${RUNTYPE_CONFIG } ${ \
                    PLANNER_CONFIG } ${MODEL_CONFIG } ${START_CONFIG } ${    \
                    TARGET_CONFIG } ${SPACE_CONFIG } ${STOP_CRITERIAS }

#PLANNER_CONFIG = "--planner-alg sba_star --knn-method bf4"
#echo_and_run./ run_CRS_planner ${OUTPUT_DEST_CONFIG } ${RUNTYPE_CONFIG } ${ \
                    PLANNER_CONFIG } ${MODEL_CONFIG } ${START_CONFIG } ${    \
                    TARGET_CONFIG } ${SPACE_CONFIG } ${STOP_CRITERIAS }

#PLANNER_CONFIG = "--planner-alg sba_star --knn-method bf4 --with-bnb"
#echo_and_run./ run_CRS_planner ${OUTPUT_DEST_CONFIG } ${RUNTYPE_CONFIG } ${ \
                    PLANNER_CONFIG } ${MODEL_CONFIG } ${START_CONFIG } ${    \
                    TARGET_CONFIG } ${SPACE_CONFIG } ${STOP_CRITERIAS }

        PLANNER_CONFIG =
            "--planner-alg sba_star --knn-method bf4 --relaxation-factor "
            "2" echo_and_run./
            run_CRS_planner ${OUTPUT_DEST_CONFIG} ${RUNTYPE_CONFIG} ${
                PLANNER_CONFIG} ${MODEL_CONFIG} ${START_CONFIG} ${
                TARGET_CONFIG} ${SPACE_CONFIG} ${STOP_CRITERIAS} mv pp_results /
            CRS_static / scene2 /
            sba_star__any_adj_list_bf4_solutions.txt pp_results / CRS_static /
            scene2 /
            sba_star__any_adj_list_bf4_r02_solutions.txt mv pp_results /
            CRS_static / scene2 /
            sba_star__any_adj_list_bf4_times.txt pp_results / CRS_static /
            scene2 /
            sba_star__any_adj_list_bf4_r02_times.txt

                PLANNER_CONFIG =
                "--planner-alg sba_star --knn-method bf4 --relaxation-factor "
                "5" echo_and_run./
                run_CRS_planner ${OUTPUT_DEST_CONFIG} ${RUNTYPE_CONFIG} ${
                    PLANNER_CONFIG} ${MODEL_CONFIG} ${START_CONFIG} ${
                    TARGET_CONFIG} ${SPACE_CONFIG} ${
                    STOP_CRITERIAS} mv pp_results /
                CRS_static / scene2 /
                sba_star__any_adj_list_bf4_solutions.txt pp_results /
                CRS_static / scene2 /
                sba_star__any_adj_list_bf4_r05_solutions.txt mv pp_results /
                CRS_static / scene2 /
                sba_star__any_adj_list_bf4_times.txt pp_results / CRS_static /
                scene2 /
                sba_star__any_adj_list_bf4_r05_times.txt

                    PLANNER_CONFIG =
                    "--planner-alg sba_star --knn-method bf4 "
                    "--relaxation-factor 10" echo_and_run./
                    run_CRS_planner ${OUTPUT_DEST_CONFIG} ${RUNTYPE_CONFIG} ${
                        PLANNER_CONFIG} ${MODEL_CONFIG} ${START_CONFIG} ${
                        TARGET_CONFIG} ${SPACE_CONFIG} ${
                        STOP_CRITERIAS} mv pp_results /
                    CRS_static / scene2 /
                    sba_star__any_adj_list_bf4_solutions.txt pp_results /
                    CRS_static / scene2 /
                    sba_star__any_adj_list_bf4_r10_solutions.txt mv pp_results /
                    CRS_static / scene2 /
                    sba_star__any_adj_list_bf4_times.txt pp_results /
                    CRS_static / scene2 /
                    sba_star__any_adj_list_bf4_r10_times.txt

#PLANNER_CONFIG = \
    "--planner-alg sba_star --knn-method bf4 --relaxation-factor 10 --with-bnb"
#echo_and_run./ run_CRS_planner ${OUTPUT_DEST_CONFIG } ${RUNTYPE_CONFIG } ${ \
                    PLANNER_CONFIG } ${MODEL_CONFIG } ${START_CONFIG } ${    \
                    TARGET_CONFIG } ${SPACE_CONFIG } ${STOP_CRITERIAS }

                        PLANNER_CONFIG =
                        "--planner-alg sba_star --knn-method bf4 "
                        "--relaxation-factor 2 --with-voronoi-pull "
                        "--sa-temperature 5" echo_and_run./
                        run_CRS_planner ${OUTPUT_DEST_CONFIG} ${
                            RUNTYPE_CONFIG} ${PLANNER_CONFIG} ${MODEL_CONFIG} ${
                            START_CONFIG} ${TARGET_CONFIG} ${SPACE_CONFIG} ${
                            STOP_CRITERIAS} mv pp_results /
                        CRS_static / scene2 /
                        sba_star__any_sa_adj_list_bf4_solutions.txt pp_results /
                        CRS_static / scene2 /
                        sba_star__any_sa_adj_list_bf4_r02_t05_solutions.txt mv
                            pp_results /
                        CRS_static / scene2 /
                        sba_star__any_sa_adj_list_bf4_times.txt pp_results /
                        CRS_static / scene2 /
                        sba_star__any_sa_adj_list_bf4_r02_t05_times.txt

                            PLANNER_CONFIG =
                            "--planner-alg sba_star --knn-method bf4 "
                            "--relaxation-factor 5 --with-voronoi-pull "
                            "--sa-temperature 5" echo_and_run./
                            run_CRS_planner ${OUTPUT_DEST_CONFIG} ${
                                RUNTYPE_CONFIG} ${PLANNER_CONFIG} ${
                                MODEL_CONFIG} ${START_CONFIG} ${
                                TARGET_CONFIG} ${SPACE_CONFIG} ${
                                STOP_CRITERIAS} mv pp_results /
                            CRS_static / scene2 /
                            sba_star__any_sa_adj_list_bf4_solutions
                                .txt pp_results /
                            CRS_static / scene2 /
                            sba_star__any_sa_adj_list_bf4_r05_t05_solutions
                                .txt mv pp_results /
                            CRS_static / scene2 /
                            sba_star__any_sa_adj_list_bf4_times.txt pp_results /
                            CRS_static / scene2 /
                            sba_star__any_sa_adj_list_bf4_r05_t05_times.txt

                                PLANNER_CONFIG =
                                "--planner-alg sba_star --knn-method bf4 "
                                "--relaxation-factor 10 --with-voronoi-pull "
                                "--sa-temperature 5" echo_and_run./
                                run_CRS_planner ${OUTPUT_DEST_CONFIG} ${
                                    RUNTYPE_CONFIG} ${PLANNER_CONFIG} ${
                                    MODEL_CONFIG} ${START_CONFIG} ${
                                    TARGET_CONFIG} ${SPACE_CONFIG} ${
                                    STOP_CRITERIAS} mv pp_results /
                                CRS_static / scene2 /
                                sba_star__any_sa_adj_list_bf4_solutions
                                    .txt pp_results /
                                CRS_static / scene2 /
                                sba_star__any_sa_adj_list_bf4_r10_t05_solutions
                                    .txt mv pp_results /
                                CRS_static / scene2 /
                                sba_star__any_sa_adj_list_bf4_times
                                    .txt pp_results /
                                CRS_static / scene2 /
                                sba_star__any_sa_adj_list_bf4_r10_t05_times.txt

#PLANNER_CONFIG = \
    "--planner-alg sba_star --knn-method bf4 --relaxation-factor 10 --with-voronoi-pull --sa-temperature 10"
#echo_and_run./ run_CRS_planner ${OUTPUT_DEST_CONFIG } ${RUNTYPE_CONFIG } ${ \
                    PLANNER_CONFIG } ${MODEL_CONFIG } ${START_CONFIG } ${    \
                    TARGET_CONFIG } ${SPACE_CONFIG } ${STOP_CRITERIAS }
#mv pp_results / CRS_static / scene2 /                                    \
    sba_star__any_sa_adj_list_bf4_solutions.txt pp_results / CRS_static / \
    scene2 / sba_star__any_sa_adj_list_bf4_r10_t10_solutions.txt
#mv pp_results / CRS_static / scene2 /                                         \
    sba_star__any_sa_adj_list_bf4_times.txt pp_results / CRS_static / scene2 / \
    sba_star__any_sa_adj_list_bf4_r10_t10_times.txt

#PLANNER_CONFIG = \
    "--planner-alg sba_star --knn-method bf4 --relaxation-factor 10 --with-voronoi-pull --sa-temperature 20"
#echo_and_run./ run_CRS_planner ${OUTPUT_DEST_CONFIG } ${RUNTYPE_CONFIG } ${ \
                    PLANNER_CONFIG } ${MODEL_CONFIG } ${START_CONFIG } ${    \
                    TARGET_CONFIG } ${SPACE_CONFIG } ${STOP_CRITERIAS }
#mv pp_results / CRS_static / scene2 /                                    \
    sba_star__any_sa_adj_list_bf4_solutions.txt pp_results / CRS_static / \
    scene2 / sba_star__any_sa_adj_list_bf4_r10_t20_solutions.txt
#mv pp_results / CRS_static / scene2 /                                         \
    sba_star__any_sa_adj_list_bf4_times.txt pp_results / CRS_static / scene2 / \
    sba_star__any_sa_adj_list_bf4_r10_t20_times.txt

#PLANNER_CONFIG = \
    "--planner-alg sba_star --knn-method bf4 --relaxation-factor 10 --with-voronoi-pull --sa-temperature 5 --with-bnb"
#echo_and_run./ run_CRS_planner ${OUTPUT_DEST_CONFIG } ${RUNTYPE_CONFIG } ${ \
                    PLANNER_CONFIG } ${MODEL_CONFIG } ${START_CONFIG } ${    \
                    TARGET_CONFIG } ${SPACE_CONFIG } ${STOP_CRITERIAS }
