
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

// This header has to be the first included because it will not work otherwise (really bad programming!!!)
#include <libplayerc++/playerc++.h>

#include "pptester_impl.h"

#include <QApplication>
#include <QFileDialog>
#include <QMessageBox>

#include "graph_alg/rrt_test_world.hpp"

PPTestWindow::PPTestWindow(QWidget* parent, Qt::WindowFlags flags)
    : QMainWindow(parent, flags),
      world_map_pixmap(nullptr),
      world_map_scene(nullptr),
      rrt_prop_diag(this),
      robot_running(false),
      robot_exec_thread(nullptr) {
  setupUi(this);

  connect(actionRun, SIGNAL(triggered()), this, SLOT(startPlanner()));
  connect(actionStop_4, SIGNAL(triggered()), this, SLOT(stopPlanner()));
  connect(actionLoad_World, SIGNAL(triggered()), this,
          SLOT(launchPlayerStage()));
  connect(actionLoad_Planner, SIGNAL(triggered()), this,
          SLOT(loadPathPlanner()));
  connect(actionLoad_Map, SIGNAL(triggered()), this, SLOT(loadWorldMap()));
  connect(actionClose_Map, SIGNAL(triggered()), this,
          SLOT(closeTestScenario()));
  connect(actionProperties, SIGNAL(triggered()), this,
          SLOT(setPlannerProperties()));
  connect(actionSave_Results, SIGNAL(triggered()), this, SLOT(saveResults()));
  connect(actionStart_Robot, SIGNAL(triggered()), this, SLOT(startRobot()));
  connect(actionStop_Robot, SIGNAL(triggered()), this, SLOT(stopRobot()));

  actionStop_4->setEnabled(false);
  // actionClose_Map->setEnabled(false);
  actionStart_Robot->setEnabled(false);
  actionStop_Robot->setEnabled(false);

  rrt_prop.setupUi(&rrt_prop_diag);

  connect(rrt_prop.buttonBox, SIGNAL(accepted()), this,
          SLOT(acceptRRTProperties()));

  world_map_scene = new QGraphicsScene(map_view);
  map_view->setScene(world_map_scene);
};

PPTestWindow::~PPTestWindow() {

  for (unsigned int i = 0; i < result_tabs_widgets.size(); ++i)
    tabWidget->removeTab(i + 1);

  for (std::vector<Ui::PPResultView*>::iterator it = result_tabs_views.begin();
       it != result_tabs_views.end(); ++it) {
    delete (*it)->canvas->scene();
    delete *it;
  };
  for (std::vector<QWidget*>::iterator it = result_tabs_widgets.begin();
       it != result_tabs_widgets.end(); ++it)
    delete *it;

  delete world_map_scene;
};

void PPTestWindow::loadWorldMap() {
  if (robot_running) {
    QMessageBox::information(
        this, tr("Error!"),
        tr("The robot is currently running, please stop it first!"));
    return;
  };

  closeTestScenario();

  QString filename = QFileDialog::getOpenFileName(
      this, tr("Open a World Map Image..."), QString(),
      tr("Image (*.bmp *.png *.jpg *.jpeg)"));
  if (!filename.isNull()) {

    world_map_cvimage = cv::imread(filename.toStdString());

    if (world_map_cvimage.empty()) {
      QMessageBox::information(
          this, tr("Error!"),
          tr("The image you have selected could not be loaded!"));
      return;
    };

    if (!world_map_pixmap) {
      world_map_pixmap = world_map_scene->addPixmap(QPixmap::fromImage(
          QImage(world_map_cvimage.ptr(), world_map_cvimage.size().width,
                 world_map_cvimage.size().height, QImage::Format_RGB888)
              .rgbSwapped()));
    } else {
      world_map_pixmap->pixmap() = QPixmap::fromImage(
          QImage(world_map_cvimage.ptr(), world_map_cvimage.size().width,
                 world_map_cvimage.size().height, QImage::Format_RGB888)
              .rgbSwapped());
    };
  };
};

void PPTestWindow::loadPathPlanner() {
  QMessageBox::information(this, tr("Not supported!"),
                           tr("Sorry. For this version of the application, the "
                              "only planner available is the "
                              "Rapidly-Exploring Random Tree (RRT) and it is "
                              "automatically loaded by default."));
};

void PPTestWindow::launchPlayerStage() {
  if (robot_running) {
    QMessageBox::information(
        this, tr("Error!"),
        tr("The robot is currently running, please stop it first!"));
    return;
  };
  if (best_path.empty()) {
    QMessageBox::information(this, tr("No Solution!"),
                             tr("The path-planner was either not run yet or "
                                "has not been able to find a solution!"));
    return;
  };

  QString filename = QFileDialog::getOpenFileName(
      this, tr("Open Player/Stage Configuration file..."), QString(),
      tr("Player/Stage Configuration File (*.cfg)"));
  if (!filename.isNull()) {
    std::string s = std::string("robot-player \"") + filename.toStdString() +
                    std::string("\" &");
    system(s.c_str());
  };

  actionStart_Robot->setEnabled(true);

  QMessageBox mb(QMessageBox::Information, tr("About to Start the Robot!"),
                 tr("The robot will start when you press OK."),
                 QMessageBox::Ok | QMessageBox::Cancel, this);
  int res = mb.exec();
  if (res & QMessageBox::Cancel)
    return;
  startRobot();
};

void PPTestWindow::startRobot() {
  actionStop_Robot->setEnabled(true);
  actionStart_Robot->setEnabled(false);
  robot_running = true;
  if (robot_exec_thread) {
    if (robot_exec_thread->joinable())
      robot_exec_thread->join();
    delete robot_exec_thread;
    robot_exec_thread = nullptr;
  };
  robot_exec_thread = new std::thread([this]() { executePath(); });
};

void PPTestWindow::stopRobot() {
  robot_running = false;
  if (robot_exec_thread) {
    if (robot_exec_thread->joinable())
      robot_exec_thread->join();
    delete robot_exec_thread;
    robot_exec_thread = nullptr;
  };
};

void PPTestWindow::startPlanner() {
  if (robot_running) {
    QMessageBox::information(
        this, tr("Error!"),
        tr("The robot is currently running, please stop it first!"));
    return;
  };
  if (world_map_cvimage.empty()) {
    QMessageBox::information(this, tr("World Map not loaded!"),
                             tr("You must first load a world map in order to "
                                "run the path-planner!"));
    return;
  };

  double pix_to_m =
      rrt_prop.GridWidth->value() / world_map_cvimage.size().width;

  rrt_test_world planner(
      world_map_cvimage, rrt_prop.MaxEdgeLength->value() / pix_to_m,
      rrt_prop.RobotRadius->value() / pix_to_m, rrt_prop.MaxNumVert->value(),
      rrt_prop.MaxNumSol->value(),
      [this](const cv::Mat& aImage, unsigned int aProgress) {
        updateWorldMap(aImage, aProgress);
      },
      [this](const cv::Mat& aImage, unsigned int aProgress,
             unsigned int aSolutionID, double aTotalDist) {
        showResultMap(aImage, aProgress, aSolutionID, aTotalDist);
      },
      rrt_prop.unidirCheck->isChecked(), rrt_prop.NNSearchDiv->value());

  rrt_test_world::pixel_coord p_start;
  p_start[0] =
      (rrt_prop.RobotStart_X->value() + rrt_prop.GridCenterX->value()) /
      pix_to_m;
  p_start[1] =
      world_map_cvimage.size().height -
      (rrt_prop.RobotStart_Y->value() + rrt_prop.GridCenterY->value()) /
          pix_to_m;
  planner.set_start_pos(p_start);
  rrt_test_world::pixel_coord p_goal;
  p_goal[0] = (rrt_prop.RobotGoal_X->value() + rrt_prop.GridCenterX->value()) /
              pix_to_m;
  p_goal[1] = world_map_cvimage.size().height -
              (rrt_prop.RobotGoal_Y->value() + rrt_prop.GridCenterY->value()) /
                  pix_to_m;
  planner.set_goal_pos(p_goal);

  planner.run();

  planner.get_best_solution(best_path);
};

void PPTestWindow::stopPlanner(){

};

void PPTestWindow::setPlannerProperties() {
  if (robot_running) {
    QMessageBox::information(
        this, tr("Error!"),
        tr("The robot is currently running, please stop it first!"));
    return;
  };
  rrt_prop_diag.exec();
};

void PPTestWindow::closeTestScenario() {
  if (robot_running)
    stopRobot();

  if (!world_map_cvimage.empty()) {
    world_map_cvimage.release();
  };

  world_map_scene->clear();
  world_map_pixmap = nullptr;

  best_path.clear();

  for (unsigned int i = 0; i < result_tabs_widgets.size(); ++i)
    tabWidget->removeTab(i + 1);

  for (std::vector<Ui::PPResultView*>::iterator it = result_tabs_views.begin();
       it != result_tabs_views.end(); ++it) {
    delete (*it)->canvas->scene();
    delete *it;
  };
  result_tabs_views.clear();
  result_names.clear();
  result_tabs_pixmaps.clear();
  for (std::vector<QWidget*>::iterator it = result_tabs_widgets.begin();
       it != result_tabs_widgets.end(); ++it)
    delete *it;
  result_tabs_widgets.clear();
};

void PPTestWindow::saveResults() {
  if (result_tabs_views.size() == 0) {
    QMessageBox::information(this, tr("No Results to Save!"),
                             tr("There are no recorded results to be saved!"));
    return;
  };
  QString dir = QFileDialog::getExistingDirectory(
      this, tr("Select a directory for the result files..."));
  if (dir.isEmpty())
    return;
  for (unsigned int i = 0; i < result_tabs_pixmaps.size(); ++i) {
    result_tabs_pixmaps[i]->pixmap().save(
        QString::fromStdString(dir.toStdString() + "/" + result_names[i]));
  };
};

void PPTestWindow::acceptRRTProperties() {
  // for the moment I don't see any invalid condition to check for.
  rrt_prop_diag.accept();
};

void PPTestWindow::updateWorldMap(const cv::Mat& aImage,
                                  unsigned int aProgress) {
  if (!world_map_pixmap) {
    world_map_pixmap = world_map_scene->addPixmap(
        QPixmap::fromImage(QImage(aImage.ptr(), aImage.size().width,
                                  aImage.size().height, QImage::Format_RGB888)
                               .rgbSwapped()));
  } else {
    world_map_pixmap->pixmap() =
        QPixmap::fromImage(QImage(aImage.ptr(), aImage.size().width,
                                  aImage.size().height, QImage::Format_RGB888)
                               .rgbSwapped());
  };
  map_view->repaint();
};

void PPTestWindow::showResultMap(const cv::Mat& aImage, unsigned int aProgress,
                                 unsigned int aSolutionID, double aTotalDist) {
  result_tabs_views.push_back(new Ui::PPResultView);
  result_tabs_widgets.push_back(new QWidget(tabWidget));
  {
    std::stringstream ss;
    ss << "Result " << aSolutionID;
    tabWidget->addTab(result_tabs_widgets.back(),
                      QString::fromStdString(ss.str()));
    result_tabs_views.back()->setupUi(result_tabs_widgets.back());
    ss << " with " << aProgress
       << " vertices, achieves a total traveled distance of " << aTotalDist;
    result_tabs_views.back()->label->setText(QString::fromStdString(ss.str()));
  };
  result_tabs_views.back()->canvas->setScene(
      new QGraphicsScene(result_tabs_views.back()->canvas));
  result_tabs_pixmaps.push_back(
      result_tabs_views.back()->canvas->scene()->addPixmap(
          QPixmap::fromImage(QImage(aImage.ptr(), aImage.size().width,
                                    aImage.size().height, QImage::Format_RGB888)
                                 .rgbSwapped())));
  {
    std::stringstream ss;
    ss << "rrt_result_v" << aProgress << "_s" << aSolutionID << "_d"
       << aTotalDist << ".png";
    result_names.push_back(ss.str());
  };
};

void PPTestWindow::executePath() {
  using namespace PlayerCc;

  PlayerClient robot("localhost");
  RangerProxy rp(&robot, 1);
  Position2dProxy pp(&robot, 0);

  double pix_to_m =
      rrt_prop.GridWidth->value() / world_map_cvimage.size().width;

  std::list<rrt_test_world::pixel_coord>::iterator pt_it = best_path.begin();

  while ((robot_running) && (!pp.GetStall()) && (pt_it != best_path.end())) {

    // read from the proxies
    robot.Read();

    double pd_x = (*pt_it)[0] * pix_to_m - rrt_prop.GridCenterX->value();
    double pd_y = (world_map_cvimage.size().height - (*pt_it)[1]) * pix_to_m -
                  rrt_prop.GridCenterY->value();

    double p_x = pp.GetXPos();
    double p_y = pp.GetYPos();

    double pd_t = atan2(pd_y - p_y, pd_x - p_x);

    player_pose2d_t goto_pose = {pd_x, pd_y, pd_t};
    pp.GoTo(goto_pose);

    if ((p_x - pd_x) * (p_x - pd_x) + (p_y - pd_y) * (p_y - pd_y) <
        2 * rrt_prop.RobotRadius->value() * rrt_prop.RobotRadius->value()) {
      pt_it++;
    };
  };

  if (pp.GetStall()) {
    std::cout << "Your robot has crashed or stalled!" << std::endl;
  } else if (pt_it == best_path.end()) {
    std::cout << "Your robot has successfully followed the path!" << std::endl;
  };

  robot_running = false;
  actionStart_Robot->setEnabled(true);
  actionStop_Robot->setEnabled(false);
};

int main(int argc, char** argv) {
  QApplication app(argc, argv);
  PPTestWindow window;
  window.show();

  return app.exec();
};
