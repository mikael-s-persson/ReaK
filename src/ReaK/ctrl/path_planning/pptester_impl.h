
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

#ifndef REAK_PPTESTER_IMPL_H
#define REAK_PPTESTER_IMPL_H

#include "ui_pptester.h"
#include "ui_ppresultview.h"
#include "ui_rrtproperties.h"

#include <vector>
#include <opencv/highgui.h>

#include <QGraphicsPixmapItem>
#include "graph_alg/rrt_test_world.hpp"
#include <boost/thread/thread.hpp>

class PPTestWindow : public QMainWindow, private Ui::PPTestWindow {
    Q_OBJECT
  
  public:
    PPTestWindow( QWidget * parent = 0, Qt::WindowFlags flags = 0 );
    ~PPTestWindow();
    
  private slots:
    
    void loadWorldMap();
    void loadPathPlanner();
    void launchPlayerStage();
    void startRobot();
    void stopRobot();
    void startPlanner();
    void stopPlanner();
    void setPlannerProperties();
    void closeTestScenario();
    void saveResults();
    void acceptRRTProperties();
    
  private:
    std::vector< Ui::PPResultView* > result_tabs_views;
    std::vector< QGraphicsPixmapItem* > result_tabs_pixmaps;
    std::vector< std::string > result_names;
    std::vector< QWidget* > result_tabs_widgets;
        
    cv::Mat world_map_cvimage;
    QGraphicsPixmapItem* world_map_pixmap;
    QGraphicsScene* world_map_scene;
    
    QDialog rrt_prop_diag;
    Ui::RRTProperties rrt_prop;
    
    std::list<rrt_test_world::pixel_coord> best_path;
    bool robot_running;
    boost::thread* robot_exec_thread;
    
    void updateWorldMap(const cv::Mat& aImage, unsigned int aProgress);
    
    void showResultMap(const cv::Mat& aImage, unsigned int aProgress, unsigned int aSolutionID, double aTotalDist);
    
    void executePath();
  
};



#endif














