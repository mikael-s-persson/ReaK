
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

#ifndef REAK_PPTESTER_PLANNER_HPP
#define REAK_PPTESTER_PLANNER_HPP

//This header has to be the first included because it will not work otherwise (really bad programming!!!)
#include <libplayerc++/playerc++.h>

#include "lin_alg/vect_alg.hpp"

#include "serialization/archiver.hpp"

#include <string>

#include <opencv/highgui.h>

class pptester_planner {
  protected:
    
    double mRobotCollisionRadius;
    
    ReaK::vect<double,2> mGridDim;
    ReaK::vect<double,2> mGridCenter;
    ReaK::vect<double,2> mStart;
    ReaK::vect<double,2> mGoal;
    
    bool mRobotRunning;
    
    std::list< ReaK::vect<double,2> > mBestPath;
    
    boost::function< void() > mCompletedMotionCallback;
    
  public:
    
    
    bool& RobotRunning() { return mRobotRunning; };
    bool RobotRunning() const { return mRobotRunning; };
    
    pptester_planner() : mRobotCollisionRadius(0.4), 
                         mGridDim(16.0,16.0),
                         mGridCenter(8.0,8.0),
                         mStart(-6.43,-5.89),
                         mGoal(7.0,7.0),
                         mRobotRunning(false),
                         mBestPath() { };
    
    virtual bool run(const cv::Mat& aWorldMap) = 0;
    
    virtual std::string getName() const = 0;
    
    void executePath(boost::function< void() > aCompletedMotionCallback) {
      mCompletedMotionCallback = aCompletedMotionCallback;
      this->executePath_impl();
    };


  protected:

    virtual void executePath_impl() {
      using namespace PlayerCc;

      PlayerClient    robot("localhost");
      RangerProxy     rp(&robot,1);
      Position2dProxy pp(&robot,0);
  
      std::list< ReaK::vect<double,2> >::iterator pt_it = mBestPath.begin();
  
      while((mRobotRunning) && (!pp.GetStall()) && (pt_it != mBestPath.end())) {
    
        // read from the proxies
        robot.Read();
    
        double pd_x = (*pt_it)[0];// * pix_to_m - rrt_prop.GridCenterX->value();
        double pd_y = (*pt_it)[1];//(world_map_cvimage.size().height - (*pt_it)[1]) * pix_to_m - rrt_prop.GridCenterY->value();
    
        double p_x = pp.GetXPos();
        double p_y = pp.GetYPos();
    
        double pd_t = atan2(pd_y - p_y, pd_x - p_x);
    
        player_pose2d_t goto_pose = {pd_x,pd_y,pd_t};
        pp.GoTo(goto_pose);
    
        if( (p_x - pd_x) * (p_x - pd_x) + (p_y - pd_y) * (p_y - pd_y) < 2*mRobotCollisionRadius*mRobotCollisionRadius) {
          pt_it++;
        };
       
      };
  
      if(pp.GetStall()) {
        std::cout << "Your robot has crashed or stalled!" << std::endl;
      } else if(pt_it == mBestPath.end()) {
        std::cout << "Your robot has successfully followed the path!" << std::endl;
      };
  
      mRobotRunning = false;
      mCompletedMotionCallback();
      //TODO: Don't forget to add this to the parent object's callback implementation.
      //actionStart_Robot->setEnabled(true);
      //actionStop_Robot->setEnabled(false);
    };
    


};




#endif











