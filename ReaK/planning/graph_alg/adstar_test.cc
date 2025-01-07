
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

#include <iostream>

#include "bagl/adjacency_list.h"
#include "bagl/astar_search.h"
#include "bagl/properties.h"

#include <functional>

#include "ReaK/math/lin_alg/vect_alg.h"
#include "ReaK/planning/graph_alg/adstar_search.h"

#include <FreeImage.h>

#include "absl/log/log.h"

//#define TESTING_ASTAR

class adstar_test_world {
 public:
  using WorldGridVertexProperties = bagl::property<
      bagl::vertex_position_t, ReaK::vect<int, 3>,
      bagl::property<
          bagl::vertex_heuristic_t, double,
          bagl::property<
              bagl::vertex_rhs_t, double,
              bagl::property<
                  bagl::vertex_key_t, ReaK::graph::adstar_key_value<double>,
                  bagl::property<
                      bagl::vertex_distance_t, double,
                      bagl::property<
                          bagl::vertex_color_t, bagl::default_color_type,
                          bagl::property<bagl::vertex_predecessor_t,
                                         bagl::adjacency_list_traits<
                                             vec_s, vec_s, bidirectional_s>::
                                             vertex_descriptor,
                                         bagl::no_property>>>>>>>;

  using WorldGridEdgeProperties =
      bagl::property<bagl::edge_weight_t, double, bagl::no_property>;

  using WorldGridType =
      bagl::adjacency_list<bagl::vec_s, bagl::vec_s, bagl::bidirectional_s,
                           WorldGridVertexProperties, WorldGridEdgeProperties>;

  using VertexType = bagl::graph_vertex_descriptor_t<WorldGridType>;
  using EdgeType = bagl::graph_edge_descriptor_t<WorldGridType>;

 private:
  FIBITMAP* world_map_image;
  FIBITMAP* world_map_output;
  int grid_width;
  int grid_height;
  int bpp;

  WorldGridType grid;
  VertexType current_pos;
  VertexType goal_pos;
  int current_time;
  int blocked_periods;
  double initial_epsilon;

  bagl::property_map_t<WorldGridType, bagl::vertex_predecessor_t> m_pred;
  bagl::property_map_t<WorldGridType, bagl::vertex_position_t> m_position;
  bagl::property_map_t<WorldGridType, bagl::vertex_distance_t> m_distance;
  bagl::property_map_t<WorldGridType, bagl::vertex_heuristic_t> m_heuristic;
  bagl::property_map_t<WorldGridType, bagl::vertex_color_t> m_color;
  bagl::property_map_t<WorldGridType, bagl::edge_weight_t> m_weight;

  double updateEdgeWeight(EdgeType e) {
    ReaK::vect<int, 3> target_pos = get(m_position, target(e, grid));
    if (((current_time < 255) && (target_pos[2] <= current_time)) ||
        ((current_time >= 255) && (target_pos[2] < 255))) {
      double old_w = get(m_weight, e);
      put(m_weight, e, (255.0 - target_pos[2]) * 10);
      return (255.0 - target_pos[2]) * 10 - old_w;
    } else {
      return 0.0;
    }
  }
  void initEdgeWeight(EdgeType e) {
    ReaK::vect<int, 3> target_pos = get(m_position, target(e, grid));
#ifdef TESTING_ASTAR
    if (target_pos[2] < 255) {
#else
    if (target_pos[2] == 0) {
#endif
      put(m_weight, e, 2550.0);
    } else {
      ReaK::vect<int, 3> source_pos = get(m_position, source(e, grid));
      if ((source_pos[0] != target_pos[0]) && (source_pos[1] != target_pos[1]))
        put(m_weight, e, std::sqrt(2.0));
      else
        put(m_weight, e, 1.0);
    }
  }

 public:
  VertexType getStartNode() {
    return goal_pos;
  }

  double adjustEpsilon(double aOldEpsilon, double aMaxWeightChange) {
    if (aMaxWeightChange > 5) {
      return initial_epsilon;
    } else {
      return (aOldEpsilon - 1.0) * 0.5 + 1.0;
    }
  }

  bool isGoalNotReached() {
    if ((current_pos == goal_pos) || (blocked_periods > 10)) {
      return false;
    } else {
      return true;
    }
  }

  template <typename EdgeIter>
  std::pair<double, EdgeIter> checkChanges(EdgeIter out_iter) {
    double max_change = 0.0;
    for (auto e : edges(grid)) {
      double ei_change = updateEdgeWeight(e);
      if (std::abs(ei_change) > 1E-3) {
        *(out_iter++) = e;
        if (std::abs(ei_change) > max_change) {
          max_change = std::abs(ei_change);
        }
      }
    }

    return std::pair<double, EdgeIter>(max_change, out_iter);
  }

  void updatePath() {

    // now update the colors of the world_image and save it.
    for (int y = 0; y < grid_height; ++y) {
      BYTE* color_bits = FreeImage_GetScanLine(world_map_output, y);
      BYTE* color_bits_orig = FreeImage_GetScanLine(world_map_image, y);

      for (int x = 0; x < grid_width; ++x) {
        VertexType current_node = vertex(y * grid_width + x, grid);
        default_color_type col = get(m_color, current_node);
        if ((color_bits_orig[FI_RGBA_RED] == color_bits_orig[FI_RGBA_GREEN]) &&
            (color_bits_orig[FI_RGBA_RED] == color_bits_orig[FI_RGBA_BLUE]) &&
            (((current_time < 255) &&
              (color_bits_orig[FI_RGBA_RED] <= current_time)) ||
             ((current_time >= 255) && (color_bits_orig[FI_RGBA_RED] < 255)))) {
          color_bits[FI_RGBA_RED] = color_bits_orig[FI_RGBA_RED];
          color_bits[FI_RGBA_GREEN] = color_bits_orig[FI_RGBA_GREEN];
          color_bits[FI_RGBA_BLUE] = color_bits_orig[FI_RGBA_BLUE];
        } else {
          if (col == white_color) {
            color_bits[FI_RGBA_RED] = 255;
            color_bits[FI_RGBA_GREEN] = 255;
            color_bits[FI_RGBA_BLUE] = 255;
          } else if (col == gray_color) {
            color_bits[FI_RGBA_RED] = 255;
            color_bits[FI_RGBA_GREEN] = 255;
            color_bits[FI_RGBA_BLUE] = 80;  // light yellow color
          } else if (col == black_color) {
            color_bits[FI_RGBA_RED] = 255;
            color_bits[FI_RGBA_GREEN] = 140;
            color_bits[FI_RGBA_BLUE] = 30;  // orange color
          } else if (col == green_color) {
            color_bits[FI_RGBA_RED] = 255;
            color_bits[FI_RGBA_GREEN] = 255;
            color_bits[FI_RGBA_BLUE] = 80;  // light yellow color
          } else {
            color_bits[FI_RGBA_RED] = 255;
            color_bits[FI_RGBA_GREEN] = 140;
            color_bits[FI_RGBA_BLUE] = 30;  // orange color
          }
          if (current_node == goal_pos) {
            color_bits[FI_RGBA_RED] = 0;
            color_bits[FI_RGBA_GREEN] = 255;
            color_bits[FI_RGBA_BLUE] = 0;  // green color
          } else if (current_node == current_pos) {
            color_bits[FI_RGBA_RED] = 0;
            color_bits[FI_RGBA_GREEN] = 0;
            color_bits[FI_RGBA_BLUE] = 255;  // blue color
          }
        }
        color_bits += bpp;
        color_bits_orig += bpp;
      }
    }

    VertexType v = current_pos;
    VertexType u = get(m_pred, v);
    std::set<VertexType> path;
    path.insert(v);
    double total_distance = 0;
    while ((u != goal_pos) && (path.insert(u).second)) {
      ReaK::vect<int, 3> p = get(m_position, u);
      BYTE* color_bits = FreeImage_GetScanLine(world_map_output, p[1]);
      // set this color to indicate the planned path.
      color_bits += bpp * p[0];
      color_bits[FI_RGBA_RED] = 255;
      color_bits[FI_RGBA_GREEN] = 0;
      color_bits[FI_RGBA_BLUE] = 0;
      total_distance += get(m_weight, edge(u, v, grid).first);
      v = u;
      u = get(m_pred, v);
    }
    total_distance += get(m_weight, edge(u, v, grid).first);
    path.clear();
    std::stringstream ss;
    ss << "test_adstar_results/" << initial_epsilon << "_" << std::setfill('0')
       << std::setw(5) << current_time << "_" << total_distance << ".bmp";
    FreeImage_Save(FIF_BMP, world_map_output, ss.str().c_str(), BMP_DEFAULT);

    // now v stores the next position, so lets move there:
    double rhs_pos = std::numeric_limits<double>::infinity();
    u = get(m_pred, current_pos);
    for (auto e : in_edges(current_pos, grid)) {
      double rhs_tmp = get(m_weight, e) + get(m_distance, source(e, grid));
      if (rhs_tmp < rhs_pos) {
        u = source(e, grid);
        rhs_pos = rhs_tmp;
      }
    }
    if (get(m_position, u)[2] == 255) {
      current_pos = u;
      blocked_periods = 0;
    } else {
      ++blocked_periods;
    }
    current_time++;
    std::cout << "\rCurrently at: " << get(m_position, u)
              << " time: " << current_time << "                ";
    std::cout.flush();

    ReaK::vect<int, 3> current_coord = get(m_position, current_pos);
    for (auto w : vertices(grid)) {
      // compute the heuristic value for each node.
      ReaK::vect<int, 3> pos = get(m_position, w);
      if (w != current_pos) {
        put(m_heuristic, w,
            std::sqrt(double(pos[0] - current_coord[0]) *
                          double(pos[0] - current_coord[0]) +
                      double(pos[1] - current_coord[1]) *
                          double(pos[1] - current_coord[1])));
      } else {
        put(m_heuristic, w, 0.0);
      }
    }

    return;
  }

  adstar_test_world(FIBITMAP* aWorldMapImage, double aInitialEpsilon)
      : world_map_image(aWorldMapImage),
        grid_width(FreeImage_GetWidth(aWorldMapImage)),
        grid_height(FreeImage_GetHeight(aWorldMapImage)),
        grid(FreeImage_GetHeight(aWorldMapImage) *
             FreeImage_GetWidth(aWorldMapImage)),
        current_time(0),
        blocked_periods(0),
        initial_epsilon(aInitialEpsilon) {
    m_pred = get(vertex_predecessor, grid);
    m_position = get(vertex_position, grid);
    m_heuristic = get(vertex_heuristic, grid);
    m_color = get(vertex_color, grid);
    m_weight = get(edge_weight, grid);
    m_distance = get(vertex_distance, grid);

    world_map_image = FreeImage_ConvertTo24Bits(world_map_image);
    world_map_output = FreeImage_Clone(world_map_image);
    FreeImage_Unload(aWorldMapImage);
    if (!world_map_image) {
      LOG(ERROR)
          << "The world image could not be converted to a 24bit Bitmap image!";
      throw int(0);
    }

    bpp = FreeImage_GetLine(world_map_image) /
          FreeImage_GetWidth(world_map_image);

    for (int y = 0; y < grid_height; ++y) {
      BYTE* color_bits = FreeImage_GetScanLine(world_map_image, y);

      for (int x = 0; x < grid_width; ++x) {
        VertexType current_node = vertex(y * grid_width + x, grid);
        ReaK::vect<int, 3> pos(x, y, 255);
        if ((color_bits[FI_RGBA_RED] == color_bits[FI_RGBA_GREEN]) &&
            (color_bits[FI_RGBA_RED] == color_bits[FI_RGBA_BLUE])) {
          // this is an obstacle of discovery time == gray value.
          if (color_bits[FI_RGBA_RED] < 254)
            pos[2] = color_bits[FI_RGBA_RED];
          else
            pos[2] = 255;
        } else if ((color_bits[FI_RGBA_RED] == 0) &&
                   (color_bits[FI_RGBA_BLUE] == 255) &&
                   (color_bits[FI_RGBA_GREEN] == 0)) {
          // this is the start position.
          pos[2] = 255;
          current_pos = current_node;
        } else if ((color_bits[FI_RGBA_RED] == 0) &&
                   (color_bits[FI_RGBA_BLUE] == 0) &&
                   (color_bits[FI_RGBA_GREEN] == 255)) {
          // this is the goal position.
          pos[2] = 255;
          goal_pos = current_node;
        }
        put(m_position, current_node, pos);

        color_bits += bpp;
      }
    }

    // std::pair<EdgeType,bool> ep = add_edge(vertex(0,grid),vertex(1,grid),grid,1.0);

    FIBITMAP* hval_image = FreeImage_ConvertToGreyscale(world_map_image);
    BYTE* hval_bits = FreeImage_GetBits(hval_image);
    ReaK::vect<int, 3> current_coord = get(m_position, current_pos);
    for (auto u : vertices(grid)) {
      // compute the heuristic value for each node.
      ReaK::vect<int, 3> pos = get(m_position, u);
      if (u != current_pos) {
        put(m_heuristic, u,
            std::sqrt(double(pos[0] - current_coord[0]) *
                          double(pos[0] - current_coord[0]) +
                      double(pos[1] - current_coord[1]) *
                          double(pos[1] - current_coord[1])));
      } else {
        put(m_heuristic, u, 0.0);
      }
      *hval_bits = int(get(m_heuristic, u)) % 256;
      hval_bits++;

      // add edges in the graph and initialize their weights.
      std::pair<EdgeType, bool> ep;
      if (pos[1] > 0) {
        ep =
            add_edge(u, vertex((pos[1] - 1) * grid_width + pos[0], grid), grid);
        if (ep.second) {
          initEdgeWeight(ep.first);
        }
        if (pos[0] > 0) {
          ep = add_edge(u, vertex((pos[1] - 1) * grid_width + pos[0] - 1, grid),
                        grid);
          if (ep.second)
            initEdgeWeight(ep.first);
        }
        if (pos[0] < grid_width - 1) {
          ep = add_edge(u, vertex((pos[1] - 1) * grid_width + pos[0] + 1, grid),
                        grid);
          if (ep.second)
            initEdgeWeight(ep.first);
        }
      }
      if (pos[1] < grid_height - 1) {
        ep =
            add_edge(u, vertex((pos[1] + 1) * grid_width + pos[0], grid), grid);
        if (ep.second) {
          initEdgeWeight(ep.first);
        }
        if (pos[0] > 0) {
          ep = add_edge(u, vertex((pos[1] + 1) * grid_width + pos[0] - 1, grid),
                        grid);
          if (ep.second) {
            initEdgeWeight(ep.first);
          }
        }
        if (pos[0] < grid_width - 1) {
          ep = add_edge(u, vertex((pos[1] + 1) * grid_width + pos[0] + 1, grid),
                        grid);
          if (ep.second) {
            initEdgeWeight(ep.first);
          }
        }
      };
      if (pos[0] > 0) {
        ep = add_edge(u, vertex(pos[1] * grid_width + pos[0] - 1, grid), grid);
        if (ep.second) {
          initEdgeWeight(ep.first);
        }
      }
      if (pos[0] < grid_width - 1) {
        ep = add_edge(u, vertex(pos[1] * grid_width + pos[0] + 1, grid), grid);
        if (ep.second) {
          initEdgeWeight(ep.first);
        }
      }
    }

    FreeImage_Save(FIF_BMP, hval_image, "test_adstar_results_hval.bmp",
                   BMP_DEFAULT);
    FreeImage_Unload(hval_image);
  }

  double getHeuristicValue(VertexType u) {
    return get(m_heuristic, u);
  }

  struct visitor {
    adstar_test_world* parent;
    explicit visitor(adstar_test_world* aParent) : parent(aParent) {}

    template <typename Vertex, typename Graph>
    void initialize_vertex(Vertex u, const Graph& g) const {
      RK_UNUSED(u);
      RK_UNUSED(g);
    }
    template <typename Vertex, typename Graph>
    void discover_vertex(Vertex u, const Graph& g) const {
      RK_UNUSED(u);
      RK_UNUSED(g);
    }
    template <typename Vertex, typename Graph>
    void inconsistent_vertex(Vertex u, const Graph& g) const {
      RK_UNUSED(u);
      RK_UNUSED(g);
    }
    template <typename Vertex, typename Graph>
    void examine_vertex(Vertex u, const Graph& g) const {
      RK_UNUSED(u);
      RK_UNUSED(g);
    }
    template <typename Edge, typename Graph>
    void examine_edge(Edge e, const Graph& g) const {
      RK_UNUSED(e);
      RK_UNUSED(g);
    }
    template <typename Edge, typename Graph>
    void edge_relaxed(Edge e, const Graph& g) const {
      RK_UNUSED(e);
      RK_UNUSED(g);
    }
    template <typename Vertex, typename Graph>
    void forget_vertex(Vertex u, const Graph& g) const {
      RK_UNUSED(u);
      RK_UNUSED(g);
    }
    template <typename Vertex, typename Graph>
    void finish_vertex(Vertex u, const Graph& g) const {
      RK_UNUSED(u);
      RK_UNUSED(g);
    }
    template <typename Vertex, typename Graph>
    void recycle_vertex(Vertex u, const Graph& g) const {
      RK_UNUSED(u);
      RK_UNUSED(g);
    }
    template <typename Graph>
    void publish_path(const Graph& g) const {
      RK_UNUSED(g);
      parent->updatePath();
    }
    bool keep_going() const { return parent->isGoalNotReached; };
    template <typename EdgeIter, typename Graph>
    std::pair<double, EdgeIter> detect_edge_change(EdgeIter ei,
                                                   const Graph& g) const {
      RK_UNUSED(g);
      return parent->checkChanges(ei);
    }
    template <typename Graph>
    double adjust_epsilon(double old_eps, double w_change,
                          const Graph& g) const {
      RK_UNUSED(g);
      return parent->adjustEpsilon(old_eps, w_change);
    }
  };

  void run() {
#ifdef TESTING_ASTAR
    current_time = 255;
    astar_search(grid, goal_pos,
                 std::bind(&adstar_test_world::getHeuristicValue, this, _1),
                 default_astar_visitor(), m_pred, m_distance,
                 get(vertex_rhs, grid), m_weight, double(0.0), m_color,
                 std::less<double>(), std::plus<double>(),
                 std::numeric_limits<double>::infinity(), double(0.0));

    updatePath();

#else
    while (goal_pos != current_pos) {
      std::cout << "Replanning from scratch..." << std::endl;

      ReaK::graph::adstar_search(grid, goal_pos, m_heuristic, visitor(this),
                                 m_pred, m_distance, get(vertex_rhs, grid),
                                 get(vertex_key, grid), m_weight, m_color,
                                 initial_epsilon);
      blocked_periods = 0;
    }
#endif
  }

  ~adstar_test_world() {
    if (world_map_image) {
      FreeImage_Unload(world_map_image);
    }
    if (world_map_output) {
      FreeImage_Unload(world_map_output);
    }
  }
};

int main(int argc, char** argv) {
  if (argc < 3) {

    LOG(ERROR)
        << "Error: Arguments to the program were incorrect!" << std::endl
        << "Usage:" << std::endl
        << "\t\t./test_adstar [world_image_filename.bmp] [initial_epsilon]"
        << std::endl;

    return 1;
  };
  std::string filename(argv[1]);
  FIBITMAP* world_image =
      FreeImage_Load(FIF_BMP, filename.c_str(), BMP_DEFAULT);

  if (!world_image) {
    LOG(ERROR) << "Error: the world image file could not be loaded! Make sure "
                  "to provide a valid BMP file!"
               << std::endl;
    return 1;
  }

  std::stringstream ss(argv[2]);
  double initial_epsilon;
  ss >> initial_epsilon;

  try {
    adstar_test_world test_world(world_image, initial_epsilon);

    test_world.run();

  } catch (...) {
    LOG(ERROR)
        << "Error: An exception was thrown during the execution of this test!"
        << std::endl;
    return 1;
  }

  return 0;
}
