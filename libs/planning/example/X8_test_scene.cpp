/**
 * \file X8_test_scene.cpp
 *
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date February 2013
 */

#include <ReaK/control/controllers/IHAQR_topology.hpp>
#include <ReaK/control/controllers/MEAQR_topology.hpp>
#include <ReaK/control/systems/quadrotor_system.hpp>
#include <ReaK/mbd/models/uav_kinematics.hpp>
#include <ReaK/topologies/spaces/se3_random_samplers.hpp>

#include <ReaK/planning/path_planning/MEAQR_rrtstar_planner.hpp>
#include <ReaK/planning/path_planning/MEAQR_sbastar_planner.hpp>

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/Qt/viewers/SoQtPlaneViewer.h>
#include <Inventor/events/SoKeyboardEvent.h>
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoEventCallback.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSwitch.h>

#include <ReaK/geometry/proximity/proxy_query_model.hpp>
#include <ReaK/mbd/coin3D/oi_scene_graph.hpp>

#include <ReaK/planning/path_planning/basic_sbmp_reporters.hpp>
#include <ReaK/planning/path_planning/frame_tracer_coin3d.hpp>

#include <ReaK/geometry/shapes/coord_arrows_3D.hpp>

#include <ReaK/mbd/kte/kte_map_chain.hpp>

#include <ReaK/core/serialization/xml_archiver.hpp>
#include <ReaK/math/optimization/optim_exceptions.hpp>

using namespace ReaK;
using namespace pp;
using namespace ctrl;

typedef IHAQR_topology<quadrotor_system::state_space_type, quadrotor_system,
                       position_only_sampler>
    X8_IHAQR_space_type;
typedef MEAQR_topology<quadrotor_system::state_space_type, quadrotor_system,
                       position_only_sampler>
    X8_MEAQR_space_type;
typedef IHAQR_topology_with_CD<quadrotor_system::state_space_type,
                               quadrotor_system, position_only_sampler>
    X8_IHAQRCD_space_type;
typedef MEAQR_topology_with_CD<quadrotor_system::state_space_type,
                               quadrotor_system, position_only_sampler>
    X8_MEAQRCD_space_type;

struct planner_params {
  std::string planner_kind;
  std::size_t max_vertices;
  std::size_t prog_interval;
  std::size_t max_results;
  std::size_t knn_method;
};

struct all_robot_info {
  std::shared_ptr<kte::UAV_kinematics> builder;
  std::shared_ptr<frame_3D<double>> X8_base_frame;
  std::shared_ptr<geom::colored_model_3D> X8_geom_model;
  std::shared_ptr<geom::proxy_query_model_3D> X8_geom_proxy;
  std::shared_ptr<kte::kte_map_chain> kin_chain;
  std::shared_ptr<frame_3D<double>> target_frame;
  SoSwitch* sw_proxy_show;
  SoSwitch* sw_motion_graph;
  SoSwitch* sw_solutions;
  vect<double, 3> start_position;
  vect<double, 3> end_position;
  std::shared_ptr<quadrotor_system> X8_sys;
  std::shared_ptr<X8_IHAQRCD_space_type> X8_IHAQR_space;
  std::shared_ptr<X8_MEAQRCD_space_type> X8_MEAQR_space;
  planner_params X8_plan_params;
};

void keyboard_press_hdl(void* userData, SoEventCallback* eventCB) {
  all_robot_info* r_info = reinterpret_cast<all_robot_info*>(userData);
  const SoEvent* event = eventCB->getEvent();

  static bool proxy_show = false;
  static bool PP_enabled = false;

  vect_n<double> j_pos = r_info->builder->getJointPositions();

  if (SO_KEY_PRESS_EVENT(event, Z)) {
    r_info->target_frame->Position[2] -= 0.01;
  } else if (SO_KEY_PRESS_EVENT(event, X)) {
    r_info->target_frame->Position[2] += 0.01;
  } else if (SO_KEY_PRESS_EVENT(event, UP_ARROW)) {
    r_info->target_frame->Position[1] -= 0.01;
  } else if (SO_KEY_PRESS_EVENT(event, DOWN_ARROW)) {
    r_info->target_frame->Position[1] += 0.01;
  } else if (SO_KEY_PRESS_EVENT(event, LEFT_ARROW)) {
    r_info->target_frame->Position[0] -= 0.01;
  } else if (SO_KEY_PRESS_EVENT(event, RIGHT_ARROW)) {
    r_info->target_frame->Position[0] += 0.01;
  } else if (SO_KEY_PRESS_EVENT(event, B)) {
    r_info->target_frame->Quat *= ReaK::axis_angle<double>(
        M_PI * 0.01, ReaK::vect<double, 3>(1.0, 0.0, 0.0));
  } else if (SO_KEY_PRESS_EVENT(event, N)) {
    r_info->target_frame->Quat *= ReaK::axis_angle<double>(
        M_PI * 0.01, ReaK::vect<double, 3>(0.0, 1.0, 0.0));
  } else if (SO_KEY_PRESS_EVENT(event, M)) {
    r_info->target_frame->Quat *= ReaK::axis_angle<double>(
        M_PI * 0.01, ReaK::vect<double, 3>(0.0, 0.0, 1.0));
  } else if (SO_KEY_PRESS_EVENT(event, L)) {
    proxy_show = !proxy_show;
    r_info->sw_proxy_show->whichChild.setValue(
        (proxy_show ? SO_SWITCH_ALL : SO_SWITCH_NONE));
  } else if (SO_KEY_PRESS_EVENT(event, P)) {
    PP_enabled = true;
  };

  r_info->kin_chain->doMotion();

  if (PP_enabled) {

    SoSeparator* mg_sep = nullptr;
    std::vector<SoSeparator*> sol_seps;

    try {
      quadrotor_system::state_space_type::point_type start_state, goal_state;
      set_position(start_state, r_info->start_position);
      set_position(goal_state, r_info->end_position);

      X8_MEAQRCD_space_type::point_type start_point(start_state);
      X8_MEAQRCD_space_type::point_type goal_point(goal_state);

      path_planning_p2p_query<X8_MEAQRCD_space_type> pp_query(
          "pp_query", r_info->X8_MEAQR_space, start_point, goal_point,
          r_info->X8_plan_params.max_results);

      any_sbmp_reporter_chain<X8_MEAQRCD_space_type> report_chain;

      //       typedef frame_tracer_3D<
      //         manip_direct_kin_map,
      //         quadrotor_system::state_space_type,
      //         MEAQR_to_state_mapper > frame_reporter_type;

      typedef frame_tracer_3D<X8_MEAQR_space_type> frame_reporter_type;

      frame_reporter_type temp_reporter(
          pp::make_any_model_applicator<X8_MEAQR_space_type>(
              manip_direct_kin_map(r_info->builder), MEAQR_to_state_mapper(),
              std::shared_ptr<quadrotor_system::state_space_type>(
                  &(r_info->X8_MEAQR_space->get_state_space()),
                  null_deleter())),
          5);
      temp_reporter.add_traced_frame(r_info->target_frame);

      report_chain.add_reporter(std::ref(temp_reporter));
      report_chain.add_reporter(print_sbmp_progress<>());

      typedef MEAQR_rrtstar_planner<quadrotor_system::state_space_type,
                                    quadrotor_system, position_only_sampler>
          X8_rrtstar_planner_type;
      X8_rrtstar_planner_type X8_planner(
          r_info->X8_MEAQR_space, r_info->X8_plan_params.max_vertices,
          r_info->X8_plan_params.prog_interval,
          ADJ_LIST_MOTION_GRAPH | r_info->X8_plan_params.knn_method,
          UNIDIRECTIONAL_PLANNING, 0.1, 0.05, 3, report_chain);

      //       typedef MEAQR_sbastar_planner<
      //         quadrotor_system::state_space_type,
      //         quadrotor_system,
      //         position_only_sampler > X8_sbastar_planner_type;
      //       X8_sbastar_planner_type X8_planner(
      //         r_info->X8_MEAQR_space,
      //         r_info->X8_plan_params.max_vertices,
      //         r_info->X8_plan_params.prog_interval,
      //         ADJ_LIST_MOTION_GRAPH | r_info->X8_plan_params.knn_method,
      //         LAZY_COLLISION_CHECKING | PLAN_WITH_VORONOI_PULL,
      //         0.1, 0.05,
      //         0.25 * r_info->X8_MEAQR_space->get_max_time_horizon() *
      //         r_info->X8_MEAQR_space->get_idle_power_cost(start_point),
      //         3, report_chain);
      //       X8_planner.set_initial_density_threshold(0.0);
      //       X8_planner.set_initial_relaxation(10.0);
      //       X8_planner.set_initial_SA_temperature(5.0);

      pp_query.reset_solution_records();
      X8_planner.solve_planning_query(pp_query);

      mg_sep = temp_reporter.get_motion_graph_tracer(r_info->target_frame)
                   .get_separator();
      mg_sep->ref();

      for (std::size_t i = 0; i < temp_reporter.get_solution_count(); ++i) {
        sol_seps.push_back(
            temp_reporter.get_solution_tracer(r_info->target_frame, i)
                .get_separator());
        sol_seps.back()->ref();
      };

    } catch (ReaK::optim::infeasible_problem& e) {
      RK_UNUSED(e);
    };
    PP_enabled = false;

    // Check the motion-graph separator and solution separators
    //  add them to the switches.
    if (mg_sep) {
      r_info->sw_motion_graph->removeAllChildren();
      r_info->sw_motion_graph->addChild(mg_sep);
      mg_sep->unref();
    };

    r_info->sw_solutions->removeAllChildren();
    if (sol_seps.size()) {
      r_info->sw_solutions->addChild(sol_seps[0]);
      sol_seps[0]->unref();
    };
  };
};

int main(int argc, char** argv) {
  using namespace ReaK;

  if (argc < 4) {
    std::cout << "Usage: ./X8_test_scene [world_filename.xml] "
                 "[quadrotor_spaces.xml] [planner_parameters.xml]"
              << std::endl
              << "\t world_filename.xml: \tThe XML filename of the world "
                 "geometry (including start and end position of "
                 "planning)."
              << std::endl
              << "\t quadrotor_spaces.xml: \tThe XML filename of the "
                 "definitions of the quad-rotor IHAQR and MEAQR topologies."
              << std::endl
              << "\t planner_parameters.xml: \tThe XML filename of the set of "
                 "parameters for the planner."
              << std::endl;
    return 1;
  };

  all_robot_info r_info;

  r_info.builder =
      std::shared_ptr<kte::UAV_kinematics>(new kte::UAV_kinematics());
  r_info.kin_chain = r_info.builder->getKTEChain();

  {
    serialization::xml_iarchive in("models/X8_quadrotor.xml");
    in >> r_info.X8_base_frame >> r_info.X8_geom_model >> r_info.X8_geom_proxy;
  };
  r_info.X8_base_frame->Parent = r_info.builder->getDependentFrame3D(0)->mFrame;

  r_info.target_frame = r_info.builder->getFrame3D(0);
  r_info.target_frame->Position = vect<double, 3>(0.0, 0.0, -1.0);

  std::shared_ptr<geom::colored_model_3D> world_geom_model;
  std::shared_ptr<geom::proxy_query_model_3D> world_geom_proxy;
  {
    std::string world_filename(argv[1]);
    serialization::xml_iarchive in(world_filename);
    in >> world_geom_model >> world_geom_proxy >> r_info.start_position >>
        r_info.end_position;
  };

  {
    std::string X8spaces_filename(argv[2]);
    serialization::xml_iarchive in(X8spaces_filename);
    std::shared_ptr<X8_IHAQR_space_type> tmp_IHAQR_space;
    std::shared_ptr<X8_MEAQR_space_type> tmp_MEAQR_space;
    in >> r_info.X8_sys >> tmp_IHAQR_space >> tmp_MEAQR_space;
    r_info.X8_IHAQR_space = std::shared_ptr<X8_IHAQRCD_space_type>(
        new X8_IHAQRCD_space_type(*tmp_IHAQR_space, r_info.builder));
    r_info.X8_MEAQR_space = std::shared_ptr<X8_MEAQRCD_space_type>(
        new X8_MEAQRCD_space_type(*tmp_MEAQR_space, r_info.builder));
  };

  std::shared_ptr<geom::proxy_query_pair_3D> X8_to_world_proxy(
      new geom::proxy_query_pair_3D("X8_to_world_proxy", r_info.X8_geom_proxy,
                                    world_geom_proxy));

  r_info.X8_IHAQR_space->m_proxy_env_3D.push_back(X8_to_world_proxy);
  r_info.X8_MEAQR_space->m_proxy_env_3D.push_back(X8_to_world_proxy);

  r_info.X8_plan_params.planner_kind = "rrt_star";
  r_info.X8_plan_params.max_vertices = 1000;
  r_info.X8_plan_params.prog_interval = 100;
  r_info.X8_plan_params.max_results = 50;
  r_info.X8_plan_params.knn_method = ReaK::pp::LINEAR_SEARCH_KNN;
  {
    std::string planner_param_filename(argv[3]);
    serialization::xml_iarchive in(planner_param_filename);
    in& RK_SERIAL_LOAD_WITH_ALIAS("planner_kind",
                                  (r_info.X8_plan_params.planner_kind)) &
        RK_SERIAL_LOAD_WITH_ALIAS("max_vertices",
                                  (r_info.X8_plan_params.max_vertices)) &
        RK_SERIAL_LOAD_WITH_ALIAS("prog_interval",
                                  (r_info.X8_plan_params.prog_interval)) &
        RK_SERIAL_LOAD_WITH_ALIAS("max_results",
                                  (r_info.X8_plan_params.max_results)) &
        RK_SERIAL_LOAD_WITH_ALIAS("knn_method",
                                  (r_info.X8_plan_params.knn_method));
  };

  {
    QWidget* mainwin = SoQt::init(argc, argv, argv[0]);

    {
      r_info.kin_chain->doMotion();

      geom::oi_scene_graph sg;

      sg.setCharacteristicLength(10.0);

      sg << (*r_info.X8_geom_model);
      sg << *world_geom_model;
      sg << geom::coord_arrows_3D(
          "base_arrows",
          std::shared_ptr<frame_3D<double>>(new frame_3D<double>()),
          pose_3D<double>(), 0.3);

      SoSeparator* root = new SoSeparator;
      root->ref();

      root->addChild(sg.getSceneGraph());

      SoEventCallback* keypressCB = new SoEventCallback;
      keypressCB->addEventCallback(SoKeyboardEvent::getClassTypeId(),
                                   keyboard_press_hdl, &r_info);
      root->addChild(keypressCB);

      r_info.sw_proxy_show = new SoSwitch();
      geom::oi_scene_graph sg_proxy;

      sg_proxy << r_info.X8_geom_proxy->getShape(0);
      r_info.sw_proxy_show->addChild(sg_proxy.getSceneGraph());
      r_info.sw_proxy_show->whichChild.setValue(SO_SWITCH_NONE);

      root->addChild(r_info.sw_proxy_show);

      r_info.sw_motion_graph = new SoSwitch();
      r_info.sw_motion_graph->whichChild.setValue(SO_SWITCH_ALL);
      root->addChild(r_info.sw_motion_graph);

      r_info.sw_solutions = new SoSwitch();
      r_info.sw_solutions->whichChild.setValue(SO_SWITCH_ALL);
      root->addChild(r_info.sw_solutions);

      sg.enableAnchorUpdates();
      sg_proxy.enableAnchorUpdates();

      // Use one of the convenient SoQt viewer classes.
      // SoQtPlaneViewer * eviewer = new SoQtPlaneViewer(mainwin);
      SoQtExaminerViewer* eviewer = new SoQtExaminerViewer(mainwin);
      eviewer->setSceneGraph(root);
      eviewer->show();

      // Pop up the main window.
      SoQt::show(mainwin);
      // Loop until exit.
      SoQt::mainLoop();

      sg_proxy.disableAnchorUpdates();
      sg.disableAnchorUpdates();

      // delete cone_rot_anim;
      // Clean up resources.
      delete eviewer;
      root->unref();
    };

    SoQt::done();
  };

  return 0;
};
