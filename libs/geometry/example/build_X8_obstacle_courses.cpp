/**
 * \file build_X8_obstacle_courses.cpp
 *
 * This application constructs the KTE-based kinematics models for the CRS A465 manipulator.
 * The inertial information is phony and is only there for completeness of the model but do
 * not reflect the actual inertial information of the CRS A465.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date February 2013
 */

#include <ReaK/geometry/shapes/box.hpp>
#include <ReaK/geometry/shapes/capped_cylinder.hpp>
#include <ReaK/geometry/shapes/colored_model.hpp>
#include <ReaK/geometry/shapes/coord_arrows_3D.hpp>
#include <ReaK/geometry/shapes/plane.hpp>

#include <ReaK/geometry/proximity/proxy_query_model.hpp>

#include <ReaK/core/serialization/xml_archiver.hpp>
#include <memory>

int main(int argc, char** argv) {

  using namespace ReaK;
  using namespace geom;

  std::shared_ptr<colored_model_3D> one_building =
      std::make_shared<colored_model_3D>("one_building_render");
  std::shared_ptr<proxy_query_model_3D> one_building_proxy =
      std::make_shared<proxy_query_model_3D>("one_building_proxy");
  {

    std::shared_ptr<plane> floor = std::make_shared<plane>(
        "floor", std::shared_ptr<pose_3D<double>>(),
        pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                        vect<double, 3>(2.5, 2.5, 0.0),
                        quaternion<double>::xrot(M_PI).getQuaternion()),
        vect<double, 2>(5.0, 5.0));

    std::shared_ptr<box> building = std::make_shared<box>(
        "building", std::shared_ptr<pose_3D<double>>(),
        pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                        vect<double, 3>(2.5, 2.5, -2.5), quaternion<double>()),
        vect<double, 3>(1.0, 1.0, 5.0));

    (*one_building)
        .addElement(color(0, 0, 0), std::make_shared<coord_arrows_3D>(
                                        "global_frame_arrows",
                                        std::shared_ptr<pose_3D<double>>(),
                                        pose_3D<double>(), 0.3))
        .addElement(color(0.3, 0.3, 0.3), floor)
        .addElement(color(0.6, 0.6, 0.0), building);

    (*one_building_proxy).addShape(floor).addShape(building);

    vect<double, 3> start_position(0.75, 0.75, -1.0);
    vect<double, 3> end_position(4.25, 4.25, -3.0);

    serialization::xml_oarchive out("models/one_building.xml");

    out << one_building << one_building_proxy << start_position << end_position;
  };

  std::shared_ptr<colored_model_3D> window_crossing =
      std::make_shared<colored_model_3D>("window_crossing_render");
  std::shared_ptr<proxy_query_model_3D> window_crossing_proxy =
      std::make_shared<proxy_query_model_3D>("window_crossing_proxy");
  {

    std::shared_ptr<plane> floor = std::make_shared<plane>(
        "floor", std::shared_ptr<pose_3D<double>>(),
        pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                        vect<double, 3>(5.0, 5.0, 0.0),
                        quaternion<double>::xrot(M_PI).getQuaternion()),
        vect<double, 2>(10.0, 10.0));

    std::shared_ptr<box> wall1 = std::make_shared<box>(
        "wall1", std::shared_ptr<pose_3D<double>>(),
        pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                        vect<double, 3>(3.0, 1.5, -5.0), quaternion<double>()),
        vect<double, 3>(0.2, 3.0, 10.0));

    std::shared_ptr<box> wall2 = std::make_shared<box>(
        "wall2", std::shared_ptr<pose_3D<double>>(),
        pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                        vect<double, 3>(3.0, 4.5, -2.5), quaternion<double>()),
        vect<double, 3>(0.2, 3.0, 5.0));

    std::shared_ptr<box> wall3 = std::make_shared<box>(
        "wall3", std::shared_ptr<pose_3D<double>>(),
        pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                        vect<double, 3>(3.0, 4.5, -9.0), quaternion<double>()),
        vect<double, 3>(0.2, 3.0, 2.0));

    std::shared_ptr<box> wall4 = std::make_shared<box>(
        "wall4", std::shared_ptr<pose_3D<double>>(),
        pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                        vect<double, 3>(3.0, 8.0, -5.0), quaternion<double>()),
        vect<double, 3>(0.2, 4.0, 10.0));

    std::shared_ptr<box> wall5 = std::make_shared<box>(
        "wall5", std::shared_ptr<pose_3D<double>>(),
        pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                        vect<double, 3>(7.0, 2.0, -5.0), quaternion<double>()),
        vect<double, 3>(0.2, 4.0, 10.0));

    std::shared_ptr<box> wall6 = std::make_shared<box>(
        "wall6", std::shared_ptr<pose_3D<double>>(),
        pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                        vect<double, 3>(7.0, 5.5, -1.0), quaternion<double>()),
        vect<double, 3>(0.2, 3.0, 2.0));

    std::shared_ptr<box> wall7 = std::make_shared<box>(
        "wall7", std::shared_ptr<pose_3D<double>>(),
        pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                        vect<double, 3>(7.0, 5.5, -7.5), quaternion<double>()),
        vect<double, 3>(0.2, 3.0, 5.0));

    std::shared_ptr<box> wall8 = std::make_shared<box>(
        "wall8", std::shared_ptr<pose_3D<double>>(),
        pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                        vect<double, 3>(7.0, 8.5, -5.0), quaternion<double>()),
        vect<double, 3>(0.2, 3.0, 5.0));

    (*window_crossing)
        .addElement(color(0, 0, 0), std::make_shared<coord_arrows_3D>(
                                        "global_frame_arrows",
                                        std::shared_ptr<pose_3D<double>>(),
                                        pose_3D<double>(), 0.3))
        .addElement(color(0.3, 0.3, 0.3), floor)
        .addElement(color(0.6, 0.6, 0.0), wall1)
        .addElement(color(0.6, 0.6, 0.0), wall2)
        .addElement(color(0.6, 0.6, 0.0), wall3)
        .addElement(color(0.6, 0.6, 0.0), wall4)
        .addElement(color(0.6, 0.6, 0.0), wall5)
        .addElement(color(0.6, 0.6, 0.0), wall6)
        .addElement(color(0.6, 0.6, 0.0), wall7)
        .addElement(color(0.6, 0.6, 0.0), wall8);

    (*window_crossing_proxy)
        .addShape(floor)
        .addShape(wall1)
        .addShape(wall2)
        .addShape(wall3)
        .addShape(wall4)
        .addShape(wall5)
        .addShape(wall6)
        .addShape(wall7)
        .addShape(wall8);

    vect<double, 3> start_position(0.75, 1.0, -1.0);
    vect<double, 3> end_position(9.0, 3.0, -7.0);

    serialization::xml_oarchive out("models/window_crossing.xml");

    out << window_crossing << window_crossing_proxy << start_position
        << end_position;
  };

  return 0;
};
