

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


#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/Qt/viewers/SoQtPlaneViewer.h>

#include <Inventor/nodes/SoSeparator.h>

#include "cylinder.hpp"
#include "box.hpp"
#include "oi_scene_graph.hpp"

#include "serialization/xml_archiver.hpp"

#include <iostream>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159
#endif

using namespace ReaK;

int main(int argc, char ** argv) {
  // Initializes SoQt library (and implicitly also the Coin and Qt
  // libraries). Returns a top-level / shell Qt window to use.
  QWidget * mainwin = SoQt::init(argc, argv, argv[0]);
  
#if 1
  pose_3D<double> a1 = pose_3D<double>(shared_ptr< pose_3D<double> >(), vect<double,3>(0.0,0.0,0.0), quaternion<double>(vect<double,4>(0.8,0.0,0.6,0.0)));
  pose_3D<double> a2 = pose_3D<double>(shared_ptr< pose_3D<double> >(), vect<double,3>(0.0,3.0,5.0), quaternion<double>(vect<double,4>(0.8,-0.6,0.0,0.0)));

  shared_ptr< geom::cylinder > cy1 = shared_ptr< geom::cylinder >(new geom::cylinder("cy1", shared_ptr< pose_3D<double> >(), a1, 5.0, 0.5));
  shared_ptr< geom::box > bx1 = shared_ptr< geom::box >(new geom::box("bx4", shared_ptr< pose_3D<double> >(), a2, vect<double,3>(4.0,2.0,2.0)));
  
  geom::colored_model_3D mdl;
  mdl.addElement(geom::color(0.5, 0.5, 0.0), cy1)
     .addElement(geom::color(0.9, 0.1, 0.1), bx1);
  
  {
    serialization::xml_oarchive out("testing_geom_serial.xml");
    out << mdl;
  };
  
#else
  
  geom::colored_model_3D mdl;
  {
    serialization::xml_iarchive in("testing_geom_serial.xml");
    in >> mdl;
  };
  
#endif
  
  {
    geom::oi_scene_graph sg;
    
    sg << mdl;
    
    // Use one of the convenient SoQt viewer classes.
    //SoQtPlaneViewer * eviewer = new SoQtPlaneViewer(mainwin);
    SoQtExaminerViewer * eviewer = new SoQtExaminerViewer(mainwin);
    eviewer->setSceneGraph(sg.getSceneGraph());
    eviewer->show();
    
    // Pop up the main window.
    SoQt::show(mainwin);
    // Loop until exit.
    SoQt::mainLoop();
    
    //delete cone_rot_anim;
    // Clean up resources.
    delete eviewer;
    
  };
  
  
  SoQt::done();
  
  return 0;
};


