
#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/Qt/viewers/SoQtPlaneViewer.h>
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoCone.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoCylinder.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoRotationXYZ.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/engines/SoTimeCounter.h>
#include <Inventor/nodes/SoRotor.h>
#include <Inventor/sensors/SoTimerSensor.h>
#include <Inventor/manips/SoTrackballManip.h>
#include <Inventor/SoOffscreenRenderer.h>

#include <iostream>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159
#endif

/*
static void incrementRotation(void *data, SoSensor*) {
  SoRotation* cone_rot = reinterpret_cast<SoRotation*>(data);
  SbRotation current_rot = cone_rot->rotation.getValue();
  current_rot *= SbRotation(SbVec3f(0,0,1), M_PI * 0.01);
  cone_rot->rotation.setValue(current_rot);
};*/


int main( int argc, char** argv ) {
  // Initializes SoQt library (and implicitly also the Coin and Qt
  // libraries). Returns a top-level / shell Qt window to use.
  QWidget* mainwin = SoQt::init( argc, argv, argv[0] );

  // Make a dead simple scene graph by using the Coin library, only
  // containing a single yellow cone under the scenegraph root.
  SoSeparator* root = new SoSeparator;
  root->ref();


  SoSeparator* cone_sep = new SoSeparator;

  SoBaseColor* col = new SoBaseColor;
  col->rgb = SbColor( 1, 1, 0 );
  cone_sep->addChild( col );

  SoTranslation* cone_pos = new SoTranslation;
  cone_pos->translation.setValue( -2.0, 0.0, 0.0 );
  cone_sep->addChild( cone_pos );

  cone_sep->addChild( new SoTrackballManip );
  /*SoRotation * cone_rot = new SoRotation;
  cone_sep->addChild(cone_rot);
  SoTimerSensor * cone_rot_anim = new SoTimerSensor(incrementRotation, cone_rot);
  cone_rot_anim->setInterval(0.02);
  cone_rot_anim->schedule();*/

  cone_sep->addChild( new SoCone );

  root->addChild( cone_sep );


  SoSeparator* cyl_sep = new SoSeparator;

  SoBaseColor* cyl_col = new SoBaseColor;
  cyl_col->rgb = SbColor( 1, 1, 0 );
  cyl_sep->addChild( cyl_col );

  SoTranslation* cyl_pos = new SoTranslation;
  cyl_pos->translation.setValue( 0.0, 0.0, 0.0 );
  cyl_sep->addChild( cyl_pos );

  //   SoRotation * cyl_rot = new SoRotation;
  //   cyl_rot->rotation.setValue(SbVec3f(1.0,0.0,0.0),M_PI/4.0);
  //   cyl_sep->addChild(cyl_rot);
  cyl_sep->addChild( new SoTrackballManip );

  SoCylinder* cyl_geom = new SoCylinder;
  cyl_geom->parts = SoCylinder::SIDES;
  cyl_sep->addChild( cyl_geom );

  SoTranslation* cyl_top_pos = new SoTranslation;
  cyl_top_pos->translation.setValue( 0.0, 1.0, 0.0 );
  cyl_sep->addChild( cyl_top_pos );

  SoSphere* cyl_top_geom = new SoSphere;
  cyl_sep->addChild( cyl_top_geom );

  SoTranslation* cyl_bottom_pos = new SoTranslation;
  cyl_bottom_pos->translation.setValue( 0.0, -2.0, 0.0 );
  cyl_sep->addChild( cyl_bottom_pos );

  SoSphere* cyl_bottom_geom = new SoSphere;
  cyl_sep->addChild( cyl_bottom_geom );


  root->addChild( cyl_sep );


  SoSeparator* cube_sep = new SoSeparator;

  SoBaseColor* cube_col = new SoBaseColor;
  cube_col->rgb = SbColor( 1, 1, 0 );
  cube_sep->addChild( cube_col );

  SoTranslation* cube_pos = new SoTranslation;
  cube_pos->translation.setValue( 2.0, 0.0, 0.0 );
  cube_sep->addChild( cube_pos );


  SoRotor* cube_rot_anim = new SoRotor;
  cube_rot_anim->rotation.setValue( SbVec3f( 1, 0, 0 ), 0.5 );
  cube_rot_anim->speed = 0.2;
  cube_sep->addChild( cube_rot_anim );

  SoCube* cube = new SoCube;
  cube->width = 1.0;
  cube->height = 2.0;
  cube->depth = 1.0;

  cube_sep->addChild( cube );

  root->addChild( cube_sep );

  // Use one of the convenient SoQt viewer classes.
  // SoQtPlaneViewer * eviewer = new SoQtPlaneViewer(mainwin);
  SoQtExaminerViewer* eviewer = new SoQtExaminerViewer( mainwin );
  eviewer->setSceneGraph( root );
  eviewer->show();

  // Pop up the main window.
  SoQt::show( mainwin );
  // Loop until exit.
  SoQt::mainLoop();

  // delete cone_rot_anim;
  // Clean up resources.
  delete eviewer;
  root->unref();
  SoQt::done();

  return 0;
};
