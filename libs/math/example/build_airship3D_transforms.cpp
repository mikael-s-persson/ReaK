
/*
 *    Copyright 2011 Sven Mikael Persson
 *
 *    THIS SOFTWARE IS DISTRIBUTED UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE v3 (GPLv3).
 *
 *    This file is part of ReaK.
 *
 *    ReaK is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    ReaK is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with ReaK (as LICENSE in the root folder).
 *    If not, see <http://www.gnu.org/licenses/>.
 */

#include <ReaK/core/serialization/xml_archiver.hpp>

#include <ReaK/math/kinetostatics/quat_alg.hpp>
#include <ReaK/math/kinetostatics/rotations_3D.hpp>

int main() {

  using namespace ReaK;

  // NOTE The rotation matrices below are in Yin Yang's convention (not mine!)

  // Airship2IMU  gives matrix such that  v_Airship = Airship2IMU * v_IMU

  std::array<double, 9> Airship2IMU_raw = {
      0.280265845357, 0.0, 0.959922421827, 0.0, -1.0, 0.0,
      0.959922421827, 0.0, -0.280265845357};

  rot_mat_3D<double> Airship2IMU(Airship2IMU_raw.data());

  quaternion<double> Airship2IMU_quat = quaternion<double>(Airship2IMU);

  quaternion<double> IMU_orientation = Airship2IMU_quat;

  vect<double, 3> IMU_location(-0.896665, 0.0, 0.25711);

  // Room2Global  gives matrix such that  v_room = Room2Global * v_gbl

  std::array<double, 9> Room2Global_raw = {
      0.75298919442, -0.65759795928, 0.02392064034,
      -0.6577626482, -0.7532241836,  -0.0012758604,
      0.018856608,   -0.0147733946,  -0.9997130464};

  rot_mat_3D<double> Room2Global(Room2Global_raw.data());

  quaternion<double> Room2Global_quat = quaternion<double>(Room2Global);

  quaternion<double> room_orientation = invert(Room2Global_quat);

  serialization::xml_oarchive out_file("models/airship3D_transforms.xml");

  out_file& RK_SERIAL_SAVE_WITH_NAME(IMU_orientation) &
      RK_SERIAL_SAVE_WITH_NAME(IMU_location) &
      RK_SERIAL_SAVE_WITH_NAME(room_orientation);

  return 0;
};
