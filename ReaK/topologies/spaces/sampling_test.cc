
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

using std::ptrdiff_t;
using std::size_t;

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "ReaK/math/lin_alg/vect_alg.h"
#include "ReaK/topologies/spaces/hyperball_topology.h"
#include "ReaK/topologies/spaces/hyperbox_topology.h"
#include "ReaK/topologies/spaces/line_topology.h"
#include "ReaK/topologies/spaces/time_poisson_topology.h"

void printResult(const cv::Mat& image, unsigned int vertex_num,
                 const std::string& dist) {
  std::stringstream ss;
  ss << "test_sampling_results/"
     << "result_" << dist << "_" << vertex_num << ".bmp";
  cv::imwrite(ss.str(), image);
};

void makeBlank(cv::Mat& image) {
  for (int i = 0; i < image.rows * image.cols; ++i) {
    image.ptr()[3 * i] = 255;
    image.ptr()[3 * i + 1] = 255;
    image.ptr()[3 * i + 2] = 255;
  };
};

int main(int argc, char** argv) {

  if (argc < 3) {

    std::cout << "Error: Arguments to the program were incorrect!" << std::endl
              << "Usage:" << std::endl
              << "\t\t./test_sampling [vertex_count] [output_interval]"
              << std::endl;

    return 1;
  };
  cv::Mat world_image(512, 512, CV_8UC3);

  if (world_image.empty()) {
    std::cout << "Error: the world image file could not be loaded! Make sure "
                 "to provide a valid BMP file!"
              << std::endl;
    return 1;
  };

  unsigned int vertex_count;
  {
    std::stringstream ss(argv[1]);
    ss >> vertex_count;
  };

  unsigned int output_interval;
  {
    std::stringstream ss(argv[2]);
    ss >> output_interval;
  };

  makeBlank(world_image);
  try {

    ReaK::pp::hyperbox_topology<ReaK::vect<double, 2>> topo(
        "hyperbox", ReaK::vect<double, 2>(0.01, 0.01),
        ReaK::vect<double, 2>(511.99, 511.99));

    for (unsigned int i = 0; i < vertex_count; ++i) {
      ReaK::vect<double, 2> v = get(ReaK::pp::random_sampler, topo)(topo);
      world_image.ptr()[3 * 512 * int(v[0]) + 3 * int(v[1])] = 0;
      world_image.ptr()[3 * 512 * int(v[0]) + 3 * int(v[1]) + 1] = 0;
      world_image.ptr()[3 * 512 * int(v[0]) + 3 * int(v[1]) + 2] = 255;
      if (i % output_interval == 0)
        printResult(world_image, i, topo.getName());
    };

    printResult(world_image, vertex_count, topo.getName());

  } catch (...) {
    std::cout
        << "Error: An exception was thrown during the execution of this test!"
        << std::endl;
    return 1;
  };

  makeBlank(world_image);
  try {

    ReaK::pp::hyperball_topology<ReaK::vect<double, 2>> topo(
        "hyperball", ReaK::vect<double, 2>(256.0, 256.0), 255.99);

    for (unsigned int i = 0; i < vertex_count; ++i) {
      ReaK::vect<double, 2> v = get(ReaK::pp::random_sampler, topo)(topo);
      world_image.ptr()[3 * 512 * int(v[0]) + 3 * int(v[1])] = 0;
      world_image.ptr()[3 * 512 * int(v[0]) + 3 * int(v[1]) + 1] = 0;
      world_image.ptr()[3 * 512 * int(v[0]) + 3 * int(v[1]) + 2] = 255;
      if (i % output_interval == 0)
        printResult(world_image, i, topo.getName());
    };

    printResult(world_image, vertex_count, topo.getName());

  } catch (...) {
    std::cout
        << "Error: An exception was thrown during the execution of this test!"
        << std::endl;
    return 1;
  };

  makeBlank(world_image);
  try {

    std::map<int, int> hist;

    ReaK::pp::time_poisson_topology topo("time_poisson", 5, 40.0);

    for (unsigned int i = 0; i < vertex_count; ++i) {
      if (hist.find(int(get(ReaK::pp::random_sampler, topo)(topo))) ==
          hist.end())
        hist[int(get(ReaK::pp::random_sampler, topo)(topo))] = 0;
      ++(hist[int(get(ReaK::pp::random_sampler, topo)(topo))]);
    };

    for (unsigned int i = 0; i < 512; ++i) {
      if (hist.find(i) == hist.end())
        hist[i] = 0;
      for (int j = 0; ((j <= hist[i]) && (j < 512)); ++j) {
        world_image.ptr()[3 * 512 * (511 - j) + 3 * i] = 0;
        world_image.ptr()[3 * 512 * (511 - j) + 3 * i + 1] = 0;
        world_image.ptr()[3 * 512 * (511 - j) + 3 * i + 2] = 255;
      };
    };

    printResult(world_image, vertex_count, topo.getName());

  } catch (...) {
    std::cout
        << "Error: An exception was thrown during the execution of this test!"
        << std::endl;
    return 1;
  };

  makeBlank(world_image);
  try {

    std::map<int, int> hist;

    ReaK::pp::line_segment_topology<double> topo("line_segment", 0, 511.99);

    for (unsigned int i = 0; i < vertex_count; ++i) {
      if (hist.find(int(get(ReaK::pp::random_sampler, topo)(topo))) ==
          hist.end())
        hist[int(get(ReaK::pp::random_sampler, topo)(topo))] = 0;
      ++(hist[int(get(ReaK::pp::random_sampler, topo)(topo))]);
    };

    for (unsigned int i = 0; i < 512; ++i) {
      if (hist.find(i) == hist.end())
        hist[i] = 0;
      for (int j = 0; ((j <= hist[i]) && (j < 512)); ++j) {
        world_image.ptr()[3 * 512 * (511 - j) + 3 * i] = 0;
        world_image.ptr()[3 * 512 * (511 - j) + 3 * i + 1] = 0;
        world_image.ptr()[3 * 512 * (511 - j) + 3 * i + 2] = 255;
      };
    };

    printResult(world_image, vertex_count, topo.getName());

  } catch (...) {
    std::cout
        << "Error: An exception was thrown during the execution of this test!"
        << std::endl;
    return 1;
  };

  return 0;
};
