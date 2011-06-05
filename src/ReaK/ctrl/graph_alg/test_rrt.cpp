
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
#include <ctime>
#include <cstdlib>

#include "rrt_test_world.hpp"


void printResult(const cv::Mat& image, unsigned int vertex_num, unsigned int solution_num, double dist) {
  std::stringstream ss;
  ss << "test_rrt_results/" << "result_" << solution_num << "_" << vertex_num << "_" << dist << ".bmp";
  cv::imwrite(ss.str(),image);
};


void printProgress(const cv::Mat& image, unsigned int vertex_num) {
  std::stringstream ss;
  ss << "test_rrt_results/" << "progress_" << vertex_num << ".bmp";
  cv::imwrite(ss.str(),image);
};




int main(int argc, char** argv) {
  std::srand(static_cast<unsigned int>(std::time(0)));

  if(argc < 7) {

    std::cout << "Error: Arguments to the program were incorrect!" << std::endl
              << "Usage:" << std::endl
              << "\t\t./test_rrt [world_image_filename.bmp] [edge_lengths] [vertex_count] [solution_count] [true|false : unidirectional] [nn_search_divider]" << std::endl;

    return 1;
  };
  std::string filename(argv[1]);
  cv::Mat world_image = cv::imread(filename);

  if(world_image.empty()) {
    std::cout << "Error: the world image file could not be loaded! Make sure to provide a valid BMP file!" << std::endl;
    return 1;
  };

  double max_edge_length;
  { std::stringstream ss(argv[2]);
    ss >> max_edge_length; };

  unsigned int max_vertex_count;
  { std::stringstream ss(argv[3]);
    ss >> max_vertex_count; };
    
  unsigned int solution_count;
  { std::stringstream ss(argv[4]);
    ss >> solution_count; };
    
  std::string is_unidir_str(argv[5]);
    
  unsigned int nn_search_divider;
  { std::stringstream ss(argv[6]);
    ss >> nn_search_divider; };

  try {
    rrt_test_world test_world(world_image,max_edge_length,0.5,max_vertex_count,solution_count,&printProgress,&printResult,is_unidir_str == "true",nn_search_divider);

    test_world.run();

  } catch (...) {
    std::cout << "Error: An exception was thrown during the execution of this test!" << std::endl;
    return 1;
  };
  

  return 0;
};















