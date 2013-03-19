
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

#include "ptrobot2D_test_world.hpp"

#ifdef REAK_HAS_OPENCV

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
//#include <opencv/cv.h>
//#include <opencv/highgui.h>

#endif

#include <string>

#include "path_planning/global_rng.hpp"

namespace ReaK {

namespace pp {


class ptrobot2D_test_world_impl {
#ifdef REAK_HAS_OPENCV
  private:
    cv::Mat world_map_image;
    mutable cv::Mat world_map_output;
    int bpp;
#endif
  public:
    int grid_width;
    int grid_height;
    ptrobot2D_test_world::point_type start_pos;
    ptrobot2D_test_world::point_type goal_pos;

    ptrobot2D_test_world_impl(const std::string& aFileName, double aRobotRadius = 1.0) {
#ifdef REAK_HAS_OPENCV
      world_map_image = cv::imread(aFileName);
      if(world_map_image.empty())
        throw std::ios_base::failure("Could not open the world map image '" + aFileName + "'! File is missing, empty or invalid!");

      world_map_output = world_map_image.clone();
      grid_width = world_map_image.size().width;
      grid_height = world_map_image.size().height;

      bpp = world_map_image.elemSize();

      for(int y = 0; y < grid_height; ++y) {
        uchar* color_bits = world_map_image.ptr(y);

        for(int x = 0; x < grid_width; ++x) {
          if( (color_bits[2] == 0) &&
              (color_bits[0] == 255) &&
              (color_bits[1] == 0) ) {
            //this is the start position.
            start_pos[0] = x;
            start_pos[1] = y;
          } else if( (color_bits[2] == 0) &&
                     (color_bits[0] == 0) &&
                     (color_bits[1] == 255) ) {
            //this is the goal position.
            goal_pos[0] = x;
            goal_pos[1] = y;
          };
          color_bits += bpp;
        };
      };

//       int iRobotRadius = int(std::fabs(aRobotRadius));
//       if(iRobotRadius > 0) {
//         cv::GaussianBlur(world_map_output,world_map_image,
//                          cv::Size(iRobotRadius * 2 + 1, iRobotRadius * 2 + 1),
//                          aRobotRadius
//                         );
//       };
#else
      grid_width = 500;
      grid_height = 500;
      start_pos = ptrobot2D_test_world::point_type(1.0,1.0);
      goal_pos  = ptrobot2D_test_world::point_type(498.0,498.0);
#endif
    };

    ptrobot2D_test_world_impl(const ptrobot2D_test_world_impl& rhs) :
#ifdef REAK_HAS_OPENCV
                              world_map_image(rhs.world_map_image.clone()),
                              world_map_output(rhs.world_map_output.clone()),
                              bpp(rhs.bpp),
#endif
                              grid_width(rhs.grid_width), grid_height(rhs.grid_height),
                              start_pos(rhs.start_pos), goal_pos(rhs.goal_pos) { };

    bool is_free(const ptrobot2D_test_world::point_type& p) const {
      if((p[0] < 0) || (p[0] >= grid_width) || (p[1] < 0) || (p[1] >= grid_height))
        return false;
#ifdef REAK_HAS_OPENCV
      const uchar* color_bits = world_map_image.ptr(int(p[1]));
      color_bits += bpp * int(p[0]);
      if( (color_bits[1] < 250) &&
          (color_bits[2] < 250) &&
          (color_bits[0] < 250) ) {
        return false;
      } else {
        return true;
      };
#else
      return true;
#endif
    };


    void reset_output() const {
#ifdef REAK_HAS_OPENCV
      world_map_output = world_map_image.clone();
#endif
    };

    void save_output(const std::string& aFilename) const {
#ifdef REAK_HAS_OPENCV
      cv::imwrite(aFilename, world_map_output);
#endif
    };

    void draw_pixel(const ptrobot2D_test_world::point_type& p, bool goal_path) const {
#ifdef REAK_HAS_OPENCV
      uchar* color_bits = world_map_output.ptr(int(p[1]));
      color_bits += bpp * int(p[0]);
      if(goal_path) {
        color_bits[2] = 255;
        color_bits[1] = 0;
        color_bits[0] = 0; //red color
      } else {
        color_bits[2] = 255;
        color_bits[1] = 140;
        color_bits[0] = 0; //orange color
      };
#endif
    };

};



bool ptrobot2D_test_world::is_free(const ptrobot2D_test_world::point_type& p) const {
  return pimpl->is_free(p);
};

void ptrobot2D_test_world::reset_output() const {
  pimpl->reset_output();
};

void ptrobot2D_test_world::save_output(const std::string& aFilename) const {
  pimpl->save_output(aFilename + ".bmp");
};

void ptrobot2D_test_world::draw_edge(const ptrobot2D_test_world::point_type& p_u, const ptrobot2D_test_world::point_type& p_v, bool goal_path) const {
  double dist = m_distance(p_v, p_u, m_space);
  if(dist < 1.0)
    return;
  double d = 0.0;
  while(d <= dist) {
    ptrobot2D_test_world::point_type p = m_space.move_position_toward(p_u, (d / dist), p_v);
    if(p[0] < 0) p[0] = 0;
    if(p[1] < 0) p[1] = 0;
    if(p[0] >= pimpl->grid_width) p[0] = pimpl->grid_width-1;
    if(p[1] >= pimpl->grid_height) p[1] = pimpl->grid_height-1;
    pimpl->draw_pixel(p,goal_path);
    d += 1.0;
  };
};


ptrobot2D_test_world::point_type ptrobot2D_test_world::random_point() const {
  point_type result;
  while(!pimpl->is_free(result = m_rand_sampler(m_space))) ; //output only free C-space points.
  return result;
};

double ptrobot2D_test_world::distance(const ptrobot2D_test_world::point_type& p1, const ptrobot2D_test_world::point_type& p2) const {
  if(m_distance(p2,move_position_toward(p1,1.0,p2), m_space) < std::numeric_limits< double >::epsilon())
    return m_distance(p1, p2, m_space); //if p2 is reachable from p1, use Euclidean distance.
  else
    return std::numeric_limits<double>::infinity(); //p2 is not reachable from p1.
};

double ptrobot2D_test_world::norm(const ptrobot2D_test_world::point_difference_type& dp) const {
  return m_distance(dp, m_space);
};

ptrobot2D_test_world::point_difference_type ptrobot2D_test_world::difference(const ptrobot2D_test_world::point_type& p1, const ptrobot2D_test_world::point_type& p2) const {
  return m_space.difference(p1,p2);
};

ptrobot2D_test_world::point_type ptrobot2D_test_world::origin() const {
  return m_space.origin();
};

ptrobot2D_test_world::point_type ptrobot2D_test_world::adjust(const ptrobot2D_test_world::point_type& p, const ptrobot2D_test_world::point_difference_type& dp) const {
  return move_position_toward(p, 1.0, m_space.adjust(p, dp));
};

ptrobot2D_test_world::point_type ptrobot2D_test_world::move_position_toward(const ptrobot2D_test_world::point_type& p1, double fraction, const ptrobot2D_test_world::point_type& p2) const {
  double dist = m_distance(p1, p2, m_space);
  if(dist * fraction > max_edge_length)
    fraction = max_edge_length / dist;
  double d = 1.0;
  while(d < dist * fraction) {
    if (!pimpl->is_free(m_space.move_position_toward(p1, (d / dist), p2))) {
      return m_space.move_position_toward(p1,((d - 1.0) / dist), p2);
    };
    d += 1.0;
  };
  if(fraction == 1.0) //these equal comparison are used for when exact end fractions are used.
    return p2;
  else if(fraction == 0.0)
    return p1;
  else
    return m_space.move_position_toward(p1, fraction, p2);
};

std::pair<ptrobot2D_test_world::point_type, bool> ptrobot2D_test_world::random_walk(const ptrobot2D_test_world::point_type& p_u) const {
  ptrobot2D_test_world::point_type p_rnd, p_v;
  unsigned int i = 0;
  do {
    p_rnd = m_rand_sampler(m_space);
    double dist = m_distance(p_u, p_rnd, m_space);
    p_v = move_position_toward(p_u, boost::uniform_01<global_rng_type&,double>(get_global_rng())() * max_edge_length / dist, p_rnd);
    ++i;
  } while((m_distance(p_u, p_v, m_space) < 1.0) && (i <= 10));
  if(i > 10) {
    //could not expand vertex u, then just generate a random C-free point.
    return std::make_pair(p_v, false);
  };
  return std::make_pair(p_v, true);
};

double ptrobot2D_test_world::bird_fly_to_goal(const ptrobot2D_test_world::point_type& p_u) const {
  return m_distance(p_u, pimpl->goal_pos, m_space);
};

double ptrobot2D_test_world::bird_fly_to_start(const ptrobot2D_test_world::point_type& p_u) const {
  return m_distance(pimpl->start_pos, p_u, m_space);
};

const ptrobot2D_test_world::point_type& ptrobot2D_test_world::get_start_pos() const {
  return pimpl->start_pos;
};

const ptrobot2D_test_world::point_type& ptrobot2D_test_world::get_goal_pos() const {
  return pimpl->goal_pos;
};

void ptrobot2D_test_world::set_start_pos(const ptrobot2D_test_world::point_type& aStart) {
  pimpl->start_pos = aStart;
};

void ptrobot2D_test_world::set_goal_pos(const ptrobot2D_test_world::point_type& aGoal) {
  pimpl->goal_pos = aGoal;
};


ptrobot2D_test_world::ptrobot2D_test_world() :
                                           pimpl(NULL),
                                           world_map_file_name(""),
                                           robot_radius(0.0),
                                           max_edge_length(0.0),
                                           m_space("ptrobot2D_space",
                                                   ptrobot2D_test_world::point_type(0,0),
                                                   ptrobot2D_test_world::point_type(0,0)),
                                           m_distance(get(distance_metric,m_space)),
                                           m_rand_sampler(get(random_sampler,m_space)) {
  setName("ptrobot2D_space_with_obstacles");
};

ptrobot2D_test_world::ptrobot2D_test_world(const std::string& aWorldMapImage,
                                           double aMaxEdgeLength,
                                           double aRobotRadius) :
                                           pimpl(new ptrobot2D_test_world_impl(aWorldMapImage, aRobotRadius)),
                                           world_map_file_name(aWorldMapImage),
                                           robot_radius(aRobotRadius),
                                           max_edge_length(aMaxEdgeLength),
                                           m_space("ptrobot2D_space",
                                                   ptrobot2D_test_world::point_type(0,0),
                                                   ptrobot2D_test_world::point_type(pimpl->grid_width,pimpl->grid_height)),
                                           m_distance(get(distance_metric,m_space)),
                                           m_rand_sampler(get(random_sampler,m_space)) {
  setName("ptrobot2D_space_with_obstacles");
};

ptrobot2D_test_world::ptrobot2D_test_world(const ptrobot2D_test_world& rhs) :
                                           pimpl(new ptrobot2D_test_world_impl(*rhs.pimpl)),
                                           world_map_file_name(rhs.world_map_file_name),
                                           robot_radius(rhs.robot_radius),
                                           max_edge_length(rhs.max_edge_length),
                                           m_space("ptrobot2D_space",
                                                   ptrobot2D_test_world::point_type(0,0),
                                                   ptrobot2D_test_world::point_type(pimpl->grid_width,pimpl->grid_height)),
                                           m_distance(get(distance_metric,m_space)),
                                           m_rand_sampler(get(random_sampler,m_space)) {
  setName("ptrobot2D_space_with_obstacles");
};

ptrobot2D_test_world& ptrobot2D_test_world::operator=(const ptrobot2D_test_world& rhs) {
  if(&rhs != this) {
    delete pimpl;
    pimpl = new ptrobot2D_test_world_impl(*rhs.pimpl);
    world_map_file_name = rhs.world_map_file_name;
    robot_radius = rhs.robot_radius;
    max_edge_length = rhs.max_edge_length;
    m_space = ptrobot2D_test_world::super_space_type("ptrobot2D_space",
                                                     ptrobot2D_test_world::point_type(0,0),
                                                     ptrobot2D_test_world::point_type(pimpl->grid_width,pimpl->grid_height));
    m_distance = get(distance_metric,m_space);
    m_rand_sampler = get(random_sampler,m_space);
  };
  return *this;
};

ptrobot2D_test_world::~ptrobot2D_test_world() {
  delete pimpl;
};

void RK_CALL ptrobot2D_test_world::save(serialization::oarchive& A, unsigned int) const {
  ReaK::named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(world_map_file_name)
    & RK_SERIAL_SAVE_WITH_NAME(robot_radius)
    & RK_SERIAL_SAVE_WITH_NAME(max_edge_length)
    & RK_SERIAL_SAVE_WITH_NAME(m_space)
    & RK_SERIAL_SAVE_WITH_NAME(m_distance)
    & RK_SERIAL_SAVE_WITH_NAME(m_rand_sampler);
};

void RK_CALL ptrobot2D_test_world::load(serialization::iarchive& A, unsigned int) {
  ReaK::named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(world_map_file_name)
    & RK_SERIAL_LOAD_WITH_NAME(robot_radius)
    & RK_SERIAL_LOAD_WITH_NAME(max_edge_length)
    & RK_SERIAL_LOAD_WITH_NAME(m_space)
    & RK_SERIAL_LOAD_WITH_NAME(m_distance)
    & RK_SERIAL_LOAD_WITH_NAME(m_rand_sampler);
  delete pimpl;
  pimpl = new ptrobot2D_test_world_impl(world_map_file_name, robot_radius);
};

};

};


