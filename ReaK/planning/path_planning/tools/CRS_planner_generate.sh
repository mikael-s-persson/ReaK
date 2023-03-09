#!/ bin / bash

./CRS_planner_run --generate_all_files models/CRS_static_base \
    --generate_xml --chaser_model_file models/CRS_A465.model.rkx \
    --target_model_file models/airship3D.model.rkx \
    --environment_models models/MD148_lab.geom.rkx \
    --planner_alg rrt_star

./CRS_planner_run --generate_chaser_target_env models/CRS_static_base_models \
    --generate_protobuf --chaser_model_file models/CRS_A465.model.rkx \
    --target_model_file models/airship3D.model.rkx \
    --environment_models models/MD148_lab.geom.rkx

./CRS_planner_run --generate_planner_options models/planners/rrtstar_base_planner \
    --generate_xml --planner_alg rrt_star --knn_method bf4

./CRS_planner_run --generate_planner_options models/planners/rrtstar_bnb_planner \
    --generate_xml --planner_alg rrt_star --knn_method bf4 --with_bnb

./CRS_planner_run --generate_planner_options models/planners/sbastar_base_planner \
    --generate_xml --planner_alg sba_star --knn_method bf4

./CRS_planner_run --generate_planner_options models/planners/sbastar_bnb_planner \
    --generate_xml --planner_alg sba_star --knn_method bf4 --with_bnb

./CRS_planner_run --generate_planner_options models/planners/sbastar_any_planner \
    --generate_xml --planner_alg sba_star --knn_method bf4 --relaxation_factor 10

./CRS_planner_run --generate_planner_options models/planners/sbastar_any_bnb_planner \
    --generate_xml --planner_alg sba_star --knn_method bf4 --relaxation_factor 10 --with_bnb

./CRS_planner_run --generate_planner_options models/planners/sbastar_any_sa_planner \
    --generate_xml --planner_alg sba_star --knn_method bf4 --relaxation_factor 10 \
    --with_voronoi_pull --sa_temperature 10

./CRS_planner_run --generate_planner_options models/planners/sbastar_any_sa_bnb_planner \
    --generate_xml --planner_alg sba_star --knn_method bf4 --relaxation_factor 10 \
    --with_voronoi_pull --sa_temperature 10 --with_bnb

./CRS_planner_run --generate_space_definition models/space_defs/planning_space_0l \
    --generate_xml --space_order 0 --interpolation_method linear --min_travel_dist 0.1 \
    --max_travel_dist 1.0 --output_space_order 0

./CRS_planner_run --generate_space_definition models/space_defs/planning_space_1l \
    --generate_xml --space_order 1 --interpolation_method linear \
    --min_travel_dist 0.1 --max_travel_dist 1.0 --output_space_order 1

./CRS_planner_run --generate_space_definition models/space_defs/planning_space_2l \
    --generate_xml --space_order 2 --interpolation_method linear \
    --min_travel_dist 0.1 --max_travel_dist 1.0 --output_space_order 2

./CRS_planner_run --generate_space_definition models/space_defs/planning_space_1c \
    --generate_xml --space_order 1 --interpolation_method cubic \
    --min_travel_dist 0.1 --max_travel_dist 1.0 --output_space_order 1

./CRS_planner_run --generate_space_definition models/space_defs/planning_space_2c \
    --generate_xml --space_order 2 --interpolation_method cubic \
    --min_travel_dist 0.1 --max_travel_dist 1.0 --output_space_order 2

./CRS_planner_run --generate_space_definition models/space_defs/planning_space_2q \
    --generate_xml --space_order 2 --interpolation_method quintic \
    --min_travel_dist 0.1 --max_travel_dist 1.0 --output_space_order 2

./CRS_planner_run --generate_space_definition models/space_defs/planning_space_1v \
    --generate_xml --space_order 1 --interpolation_method svp \
    --min_travel_dist 0.1 --max_travel_dist 1.0 --output_space_order 1

./CRS_planner_run --generate_space_definition models/space_defs/planning_space_2v \
    --generate_xml --space_order 2 --interpolation_method svp \
    --min_travel_dist 0.1 --max_travel_dist 1.0 --output_space_order 2

./CRS_planner_run --generate_space_definition models/space_defs/planning_space_2a \
    --generate_xml --space_order 2 --interpolation_method sap \
    --min_travel_dist 0.1 --max_travel_dist 1.0 --output_space_order 2
