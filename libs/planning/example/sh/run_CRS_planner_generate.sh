#!/bin/bash

./run_CRS_planner --generate-all-files models/CRS_static_base --generate-xml --chaser-model-file models/CRS_A465.model.rkx --target-model-file models/airship3D.model.rkx --environment-models models/MD148_lab.geom.rkx --planner-alg rrt_star

./run_CRS_planner --generate-chaser-target-env models/CRS_static_base_models --generate-protobuf --chaser-model-file models/CRS_A465.model.rkx --target-model-file models/airship3D.model.rkx --environment-models models/MD148_lab.geom.rkx

./run_CRS_planner --generate-planner-options models/planners/rrtstar_base_planner --generate-xml --planner-alg rrt_star --knn-method bf4
./run_CRS_planner --generate-planner-options models/planners/rrtstar_bnb_planner --generate-xml --planner-alg rrt_star --knn-method bf4 --with-bnb

./run_CRS_planner --generate-planner-options models/planners/sbastar_base_planner --generate-xml --planner-alg sba_star --knn-method bf4
./run_CRS_planner --generate-planner-options models/planners/sbastar_bnb_planner --generate-xml --planner-alg sba_star --knn-method bf4 --with-bnb

./run_CRS_planner --generate-planner-options models/planners/sbastar_any_planner --generate-xml --planner-alg sba_star --knn-method bf4 --relaxation-factor 10
./run_CRS_planner --generate-planner-options models/planners/sbastar_any_bnb_planner --generate-xml --planner-alg sba_star --knn-method bf4 --relaxation-factor 10 --with-bnb

./run_CRS_planner --generate-planner-options models/planners/sbastar_any_sa_planner --generate-xml --planner-alg sba_star --knn-method bf4 --relaxation-factor 10 --with-voronoi-pull --sa-temperature 10
./run_CRS_planner --generate-planner-options models/planners/sbastar_any_sa_bnb_planner --generate-xml --planner-alg sba_star --knn-method bf4 --relaxation-factor 10 --with-voronoi-pull --sa-temperature 10 --with-bnb

./run_CRS_planner --generate-space-definition models/space_defs/planning_space_0l --generate-xml --space-order 0 --interpolation-method linear --min-travel-dist 0.1 --max-travel-dist 1.0 --output-space-order 0
./run_CRS_planner --generate-space-definition models/space_defs/planning_space_1l --generate-xml --space-order 1 --interpolation-method linear --min-travel-dist 0.1 --max-travel-dist 1.0 --output-space-order 1
./run_CRS_planner --generate-space-definition models/space_defs/planning_space_2l --generate-xml --space-order 2 --interpolation-method linear --min-travel-dist 0.1 --max-travel-dist 1.0 --output-space-order 2

./run_CRS_planner --generate-space-definition models/space_defs/planning_space_1c --generate-xml --space-order 1 --interpolation-method cubic --min-travel-dist 0.1 --max-travel-dist 1.0 --output-space-order 1
./run_CRS_planner --generate-space-definition models/space_defs/planning_space_2c --generate-xml --space-order 2 --interpolation-method cubic --min-travel-dist 0.1 --max-travel-dist 1.0 --output-space-order 2

./run_CRS_planner --generate-space-definition models/space_defs/planning_space_2q --generate-xml --space-order 2 --interpolation-method quintic --min-travel-dist 0.1 --max-travel-dist 1.0 --output-space-order 2

./run_CRS_planner --generate-space-definition models/space_defs/planning_space_1v --generate-xml --space-order 1 --interpolation-method svp --min-travel-dist 0.1 --max-travel-dist 1.0 --output-space-order 1
./run_CRS_planner --generate-space-definition models/space_defs/planning_space_2v --generate-xml --space-order 2 --interpolation-method svp --min-travel-dist 0.1 --max-travel-dist 1.0 --output-space-order 2

./run_CRS_planner --generate-space-definition models/space_defs/planning_space_2a --generate-xml --space-order 2 --interpolation-method sap --min-travel-dist 0.1 --max-travel-dist 1.0 --output-space-order 2


