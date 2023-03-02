
/*
 *    Copyright 2012 Sven Mikael Persson
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

#include <ReaK/mbd/models/manip_kinematics_helper.hpp>

namespace ReaK::kte {

void manip_kin_mdl_joint_io::getJointPositions(double* result) const {

  for (std::size_t i = 0; i < model->getCoordsCount(); ++i) {
    (*result++) = model->getCoord(i)->q;
  }

  for (std::size_t i = 0; i < model->getFrames2DCount(); ++i) {
    std::shared_ptr<frame_2D<double>> i_ptr = model->getFrame2D(i);
    (*result++) = i_ptr->Position[0];
    (*result++) = i_ptr->Position[1];
    (*result++) = i_ptr->Rotation[0];
    (*result++) = i_ptr->Rotation[1];
  }

  for (std::size_t i = 0; i < model->getFrames3DCount(); ++i) {
    std::shared_ptr<frame_3D<double>> i_ptr = model->getFrame3D(i);
    (*result++) = i_ptr->Position[0];
    (*result++) = i_ptr->Position[1];
    (*result++) = i_ptr->Position[2];
    (*result++) = i_ptr->Quat[0];
    (*result++) = i_ptr->Quat[1];
    (*result++) = i_ptr->Quat[2];
    (*result++) = i_ptr->Quat[3];
  }
}

void manip_kin_mdl_joint_io::setJointPositions(const double* aJointPositions) {

  for (std::size_t i = 0; i < model->getCoordsCount(); ++i) {
    model->getCoord(i)->q = (*aJointPositions++);
  }

  for (std::size_t i = 0; i < model->getFrames2DCount(); ++i) {
    std::shared_ptr<frame_2D<double>> i_ptr = model->getFrame2D(i);
    i_ptr->Position[0] = (*aJointPositions++);
    i_ptr->Position[1] = (*aJointPositions++);
    i_ptr->Rotation = rot_mat_2D<double>(
        vect<double, 2>(aJointPositions[0], aJointPositions[1]));
    aJointPositions += 2;
  }

  for (std::size_t i = 0; i < model->getFrames3DCount(); ++i) {
    std::shared_ptr<frame_3D<double>> i_ptr = model->getFrame3D(i);
    i_ptr->Position[0] = (*aJointPositions++);
    i_ptr->Position[1] = (*aJointPositions++);
    i_ptr->Position[2] = (*aJointPositions++);
    i_ptr->Quat = quaternion<double>(
        vect<double, 4>(aJointPositions[0], aJointPositions[1],
                        aJointPositions[2], aJointPositions[3]));
    aJointPositions += 4;
  }
}

void manip_kin_mdl_joint_io::getJointVelocities(double* result) const {

  for (std::size_t i = 0; i < model->getCoordsCount(); ++i) {
    (*result++) = model->getCoord(i)->q_dot;
  }

  for (std::size_t i = 0; i < model->getFrames2DCount(); ++i) {
    std::shared_ptr<frame_2D<double>> i_ptr = model->getFrame2D(i);
    (*result++) = i_ptr->Velocity[0];
    (*result++) = i_ptr->Velocity[1];
    (*result++) = i_ptr->AngVelocity;
  }

  for (std::size_t i = 0; i < model->getFrames3DCount(); ++i) {
    std::shared_ptr<frame_3D<double>> i_ptr = model->getFrame3D(i);
    (*result++) = i_ptr->Velocity[0];
    (*result++) = i_ptr->Velocity[1];
    (*result++) = i_ptr->Velocity[2];
    (*result++) = i_ptr->AngVelocity[0];
    (*result++) = i_ptr->AngVelocity[1];
    (*result++) = i_ptr->AngVelocity[2];
  }
}

void manip_kin_mdl_joint_io::setJointVelocities(
    const double* aJointVelocities) {

  for (std::size_t i = 0; i < model->getCoordsCount(); ++i) {
    model->getCoord(i)->q_dot = (*aJointVelocities++);
  }

  for (std::size_t i = 0; i < model->getFrames2DCount(); ++i) {
    std::shared_ptr<frame_2D<double>> i_ptr = model->getFrame2D(i);
    i_ptr->Velocity[0] = (*aJointVelocities++);
    i_ptr->Velocity[1] = (*aJointVelocities++);
    i_ptr->AngVelocity = (*aJointVelocities++);
  }

  for (std::size_t i = 0; i < model->getFrames3DCount(); ++i) {
    std::shared_ptr<frame_3D<double>> i_ptr = model->getFrame3D(i);
    i_ptr->Velocity[0] = (*aJointVelocities++);
    i_ptr->Velocity[1] = (*aJointVelocities++);
    i_ptr->Velocity[2] = (*aJointVelocities++);
    i_ptr->AngVelocity[0] = (*aJointVelocities++);
    i_ptr->AngVelocity[1] = (*aJointVelocities++);
    i_ptr->AngVelocity[2] = (*aJointVelocities++);
  }
}

void manip_kin_mdl_joint_io::getJointAccelerations(double* result) const {

  for (std::size_t i = 0; i < model->getCoordsCount(); ++i) {
    (*result++) = model->getCoord(i)->q_ddot;
  }

  for (std::size_t i = 0; i < model->getFrames2DCount(); ++i) {
    std::shared_ptr<frame_2D<double>> i_ptr = model->getFrame2D(i);
    (*result++) = i_ptr->Acceleration[0];
    (*result++) = i_ptr->Acceleration[1];
    (*result++) = i_ptr->AngAcceleration;
  }

  for (std::size_t i = 0; i < model->getFrames3DCount(); ++i) {
    std::shared_ptr<frame_3D<double>> i_ptr = model->getFrame3D(i);
    (*result++) = i_ptr->Acceleration[0];
    (*result++) = i_ptr->Acceleration[1];
    (*result++) = i_ptr->Acceleration[2];
    (*result++) = i_ptr->AngAcceleration[0];
    (*result++) = i_ptr->AngAcceleration[1];
    (*result++) = i_ptr->AngAcceleration[2];
  }
}

void manip_kin_mdl_joint_io::setJointAccelerations(
    const double* aJointAccelerations) {

  for (std::size_t i = 0; i < model->getCoordsCount(); ++i) {
    model->getCoord(i)->q_ddot = (*aJointAccelerations++);
  }

  for (std::size_t i = 0; i < model->getFrames2DCount(); ++i) {
    std::shared_ptr<frame_2D<double>> i_ptr = model->getFrame2D(i);
    i_ptr->Acceleration[0] = (*aJointAccelerations++);
    i_ptr->Acceleration[1] = (*aJointAccelerations++);
    i_ptr->AngAcceleration = (*aJointAccelerations++);
  }

  for (std::size_t i = 0; i < model->getFrames3DCount(); ++i) {
    std::shared_ptr<frame_3D<double>> i_ptr = model->getFrame3D(i);
    i_ptr->Acceleration[0] = (*aJointAccelerations++);
    i_ptr->Acceleration[1] = (*aJointAccelerations++);
    i_ptr->Acceleration[2] = (*aJointAccelerations++);
    i_ptr->AngAcceleration[0] = (*aJointAccelerations++);
    i_ptr->AngAcceleration[1] = (*aJointAccelerations++);
    i_ptr->AngAcceleration[2] = (*aJointAccelerations++);
  }
}

void manip_kin_mdl_joint_io::getDependentPositions(double* result) const {

  for (std::size_t i = 0; i < model->getDependentCoordsCount(); ++i) {
    (*result++) = model->getDependentCoord(i)->mFrame->q;
  }

  for (std::size_t i = 0; i < model->getDependentFrames2DCount(); ++i) {
    std::shared_ptr<joint_dependent_frame_2D> i_ptr =
        model->getDependentFrame2D(i);
    (*result++) = i_ptr->mFrame->Position[0];
    (*result++) = i_ptr->mFrame->Position[1];
    (*result++) = i_ptr->mFrame->Rotation[0];
    (*result++) = i_ptr->mFrame->Rotation[1];
  }

  for (std::size_t i = 0; i < model->getDependentFrames3DCount(); ++i) {
    std::shared_ptr<joint_dependent_frame_3D> i_ptr =
        model->getDependentFrame3D(i);
    (*result++) = i_ptr->mFrame->Position[0];
    (*result++) = i_ptr->mFrame->Position[1];
    (*result++) = i_ptr->mFrame->Position[2];
    (*result++) = i_ptr->mFrame->Quat[0];
    (*result++) = i_ptr->mFrame->Quat[1];
    (*result++) = i_ptr->mFrame->Quat[2];
    (*result++) = i_ptr->mFrame->Quat[3];
  }
}

void manip_kin_mdl_joint_io::getDependentVelocities(double* result) const {

  for (std::size_t i = 0; i < model->getDependentCoordsCount(); ++i) {
    (*result++) = model->getDependentCoord(i)->mFrame->q_dot;
  }

  for (std::size_t i = 0; i < model->getDependentFrames2DCount(); ++i) {
    std::shared_ptr<joint_dependent_frame_2D> i_ptr =
        model->getDependentFrame2D(i);
    (*result++) = i_ptr->mFrame->Velocity[0];
    (*result++) = i_ptr->mFrame->Velocity[1];
    (*result++) = i_ptr->mFrame->AngVelocity;
  }

  for (std::size_t i = 0; i < model->getDependentFrames3DCount(); ++i) {
    std::shared_ptr<joint_dependent_frame_3D> i_ptr =
        model->getDependentFrame3D(i);
    (*result++) = i_ptr->mFrame->Velocity[0];
    (*result++) = i_ptr->mFrame->Velocity[1];
    (*result++) = i_ptr->mFrame->Velocity[2];
    (*result++) = i_ptr->mFrame->AngVelocity[0];
    (*result++) = i_ptr->mFrame->AngVelocity[1];
    (*result++) = i_ptr->mFrame->AngVelocity[2];
  }
}

void manip_kin_mdl_joint_io::getDependentAccelerations(double* result) const {

  for (std::size_t i = 0; i < model->getDependentCoordsCount(); ++i) {
    (*result++) = model->getDependentCoord(i)->mFrame->q_ddot;
  }

  for (std::size_t i = 0; i < model->getDependentFrames2DCount(); ++i) {
    std::shared_ptr<joint_dependent_frame_2D> i_ptr =
        model->getDependentFrame2D(i);
    (*result++) = i_ptr->mFrame->Acceleration[0];
    (*result++) = i_ptr->mFrame->Acceleration[1];
    (*result++) = i_ptr->mFrame->AngAcceleration;
  }

  for (std::size_t i = 0; i < model->getDependentFrames3DCount(); ++i) {
    std::shared_ptr<joint_dependent_frame_3D> i_ptr =
        model->getDependentFrame3D(i);
    (*result++) = i_ptr->mFrame->Acceleration[0];
    (*result++) = i_ptr->mFrame->Acceleration[1];
    (*result++) = i_ptr->mFrame->Acceleration[2];
    (*result++) = i_ptr->mFrame->AngAcceleration[0];
    (*result++) = i_ptr->mFrame->AngAcceleration[1];
    (*result++) = i_ptr->mFrame->AngAcceleration[2];
  }
}

void manip_kin_mdl_jac_calculator::getJacobianMatrixAndDerivativeImpl(
    mat<double, mat_structure::rectangular>& Jac,
    mat<double, mat_structure::rectangular>* JacDot) const {
  using SubMat = mat_sub_block<mat<double, mat_structure::rectangular>>;

  std::size_t m = model->getDependentVelocitiesCount();
  std::size_t n = model->getJointVelocitiesCount();
  Jac = mat<double, mat_structure::nil>(m, n);
  if (JacDot != nullptr) {
    *JacDot = mat<double, mat_structure::nil>(m, n);
  }

  std::size_t RowInd = 0;

  /****************************************************************************************
   *                             Gen Coords
   * *************************************************************************************/

  for (std::size_t i = 0; i < model->getCoordsCount(); ++i) {
    RowInd = 0;
    std::shared_ptr<gen_coord<double>> i_joint = model->getCoord(i);

    for (std::size_t j = 0; j < model->getDependentCoordsCount(); ++j) {
      std::shared_ptr<joint_dependent_gen_coord> j_dep =
          model->getDependentCoord(j);
      if (j_dep->mUpStreamJoints.find(i_joint) !=
          j_dep->mUpStreamJoints.end()) {
        SubMat subJac = sub(Jac)(range(RowInd, RowInd + 1), range(i, i + 1));
        if (JacDot != nullptr) {
          SubMat subJacDot =
              sub(*JacDot)(range(RowInd, RowInd + 1), range(i, i + 1));
          j_dep->mUpStreamJoints[i_joint]->write_to_matrices(subJac, subJacDot);
        } else {
          j_dep->mUpStreamJoints[i_joint]->write_to_matrices(subJac);
        }
      }
      RowInd++;
    }

    for (std::size_t j = 0; j < model->getDependentFrames2DCount(); ++j) {
      std::shared_ptr<joint_dependent_frame_2D> j_dep =
          model->getDependentFrame2D(j);
      if (j_dep->mUpStreamJoints.find(i_joint) !=
          j_dep->mUpStreamJoints.end()) {
        SubMat subJac = sub(Jac)(range(RowInd, RowInd + 3), range(i, i + 1));
        if (JacDot != nullptr) {
          SubMat subJacDot =
              sub(*JacDot)(range(RowInd, RowInd + 3), range(i, i + 1));
          j_dep->mUpStreamJoints[i_joint]
              ->get_jac_relative_to(j_dep->mFrame)
              .write_to_matrices(subJac, subJacDot);
        } else {
          j_dep->mUpStreamJoints[i_joint]
              ->get_jac_relative_to(j_dep->mFrame)
              .write_to_matrices(subJac);
        }
      }
      RowInd += 3;
    }

    for (std::size_t j = 0; j < model->getDependentFrames3DCount(); ++j) {
      std::shared_ptr<joint_dependent_frame_3D> j_dep =
          model->getDependentFrame3D(j);
      if (j_dep->mUpStreamJoints.find(i_joint) !=
          j_dep->mUpStreamJoints.end()) {
        SubMat subJac = sub(Jac)(range(RowInd, RowInd + 6), range(i, i + 1));
        if (JacDot != nullptr) {
          SubMat subJacDot =
              sub(*JacDot)(range(RowInd, RowInd + 6), range(i, i + 1));
          j_dep->mUpStreamJoints[i_joint]
              ->get_jac_relative_to(j_dep->mFrame)
              .write_to_matrices(subJac, subJacDot);
        } else {
          j_dep->mUpStreamJoints[i_joint]
              ->get_jac_relative_to(j_dep->mFrame)
              .write_to_matrices(subJac);
        }
      }
      RowInd += 6;
    }
  }

  /****************************************************************************************
   *                             2D Frames
   * *************************************************************************************/

  std::size_t base_i = model->getCoordsCount();
  for (std::size_t i = 0; i < model->getFrames2DCount(); ++i) {
    RowInd = 0;
    std::shared_ptr<frame_2D<double>> i_joint = model->getFrame2D(i);

    for (std::size_t j = 0; j < model->getDependentCoordsCount(); ++j) {
      std::shared_ptr<joint_dependent_gen_coord> j_dep =
          model->getDependentCoord(j);
      if (j_dep->mUpStream2DJoints.find(i_joint) !=
          j_dep->mUpStream2DJoints.end()) {
        SubMat subJac = sub(Jac)(range(RowInd, RowInd + 1),
                                 range(3 * i + base_i, 3 * i + base_i + 3));
        if (JacDot != nullptr) {
          SubMat subJacDot =
              sub(*JacDot)(range(RowInd, RowInd + 1),
                           range(3 * i + base_i, 3 * i + base_i + 3));
          j_dep->mUpStream2DJoints[i_joint]->write_to_matrices(subJac,
                                                               subJacDot);
        } else {
          j_dep->mUpStream2DJoints[i_joint]->write_to_matrices(subJac);
        }
      }
      RowInd++;
    }

    for (std::size_t j = 0; j < model->getDependentFrames2DCount(); ++j) {
      std::shared_ptr<joint_dependent_frame_2D> j_dep =
          model->getDependentFrame2D(j);
      if (j_dep->mUpStream2DJoints.find(i_joint) !=
          j_dep->mUpStream2DJoints.end()) {
        SubMat subJac = sub(Jac)(range(RowInd, RowInd + 3),
                                 range(3 * i + base_i, 3 * i + base_i + 3));
        if (JacDot != nullptr) {
          SubMat subJacDot =
              sub(*JacDot)(range(RowInd, RowInd + 3),
                           range(3 * i + base_i, 3 * i + base_i + 3));
          j_dep->mUpStream2DJoints[i_joint]
              ->get_jac_relative_to(j_dep->mFrame)
              .write_to_matrices(subJac, subJacDot);
        } else {
          j_dep->mUpStream2DJoints[i_joint]
              ->get_jac_relative_to(j_dep->mFrame)
              .write_to_matrices(subJac);
        }
      }
      RowInd += 3;
    }

    for (std::size_t j = 0; j < model->getDependentFrames3DCount(); ++j) {
      std::shared_ptr<joint_dependent_frame_3D> j_dep =
          model->getDependentFrame3D(j);
      if (j_dep->mUpStream2DJoints.find(i_joint) !=
          j_dep->mUpStream2DJoints.end()) {
        SubMat subJac = sub(Jac)(range(RowInd, RowInd + 6),
                                 range(3 * i + base_i, 3 * i + base_i + 3));
        if (JacDot != nullptr) {
          SubMat subJacDot =
              sub(*JacDot)(range(RowInd, RowInd + 6),
                           range(3 * i + base_i, 3 * i + base_i + 3));
          j_dep->mUpStream2DJoints[i_joint]
              ->get_jac_relative_to(j_dep->mFrame)
              .write_to_matrices(subJac, subJacDot);
        } else {
          j_dep->mUpStream2DJoints[i_joint]
              ->get_jac_relative_to(j_dep->mFrame)
              .write_to_matrices(subJac);
        }
      }
      RowInd += 6;
    }
  }

  /****************************************************************************************
   *                             3D Frames
   * *************************************************************************************/

  base_i = model->getCoordsCount() + 3 * model->getFrames2DCount();
  for (std::size_t i = 0; i < model->getFrames3DCount(); ++i) {
    RowInd = 0;
    std::shared_ptr<frame_3D<double>> i_joint = model->getFrame3D(i);

    for (std::size_t j = 0; j < model->getDependentCoordsCount(); ++j) {
      std::shared_ptr<joint_dependent_gen_coord> j_dep =
          model->getDependentCoord(j);
      if (j_dep->mUpStream3DJoints.find(i_joint) !=
          j_dep->mUpStream3DJoints.end()) {
        SubMat subJac = sub(Jac)(range(RowInd, RowInd + 1),
                                 range(6 * i + base_i, 6 * i + base_i + 6));
        if (JacDot != nullptr) {
          SubMat subJacDot =
              sub(*JacDot)(range(RowInd, RowInd + 1),
                           range(6 * i + base_i, 6 * i + base_i + 6));
          j_dep->mUpStream3DJoints[i_joint]->write_to_matrices(subJac,
                                                               subJacDot);
        } else {
          j_dep->mUpStream3DJoints[i_joint]->write_to_matrices(subJac);
        }
      }
      RowInd++;
    }

    for (std::size_t j = 0; j < model->getDependentFrames2DCount(); ++j) {
      std::shared_ptr<joint_dependent_frame_2D> j_dep =
          model->getDependentFrame2D(j);
      if (j_dep->mUpStream3DJoints.find(i_joint) !=
          j_dep->mUpStream3DJoints.end()) {
        SubMat subJac = sub(Jac)(range(RowInd, RowInd + 3),
                                 range(6 * i + base_i, 6 * i + base_i + 6));
        if (JacDot != nullptr) {
          SubMat subJacDot =
              sub(*JacDot)(range(RowInd, RowInd + 3),
                           range(6 * i + base_i, 6 * i + base_i + 6));
          j_dep->mUpStream3DJoints[i_joint]
              ->get_jac_relative_to(j_dep->mFrame)
              .write_to_matrices(subJac, subJacDot);
        } else {
          j_dep->mUpStream3DJoints[i_joint]
              ->get_jac_relative_to(j_dep->mFrame)
              .write_to_matrices(subJac);
        }
      }
      RowInd += 3;
    }

    for (std::size_t j = 0; j < model->getDependentFrames3DCount(); ++j) {
      std::shared_ptr<joint_dependent_frame_3D> j_dep =
          model->getDependentFrame3D(j);
      if (j_dep->mUpStream3DJoints.find(i_joint) !=
          j_dep->mUpStream3DJoints.end()) {
        SubMat subJac = sub(Jac)(range(RowInd, RowInd + 6),
                                 range(6 * i + base_i, 6 * i + base_i + 6));
        if (JacDot != nullptr) {
          SubMat subJacDot =
              sub(*JacDot)(range(RowInd, RowInd + 6),
                           range(6 * i + base_i, 6 * i + base_i + 6));
          j_dep->mUpStream3DJoints[i_joint]
              ->get_jac_relative_to(j_dep->mFrame)
              .write_to_matrices(subJac, subJacDot);
        } else {
          j_dep->mUpStream3DJoints[i_joint]
              ->get_jac_relative_to(j_dep->mFrame)
              .write_to_matrices(subJac);
        }
      }
      RowInd += 6;
    }
  }
}
}  // namespace ReaK::kte
