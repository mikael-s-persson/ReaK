
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

#include <ReaK/core/recorders/bin_recorder.hpp>
#include <ReaK/core/recorders/network_recorder.hpp>
#include <ReaK/core/recorders/ascii_recorder.hpp>
#include <ReaK/core/recorders/vector_recorder.hpp>

#include <sstream>

#include <ReaK/core/base/chrono_incl.hpp>
#include <ReaK/core/base/thread_incl.hpp>

#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE recorders
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


BOOST_AUTO_TEST_CASE( ascii_space_record_extract_test )
{
  using namespace ReaK;
  using namespace recorder;
  
  {
    std::stringstream ss;
    {
      ascii_recorder output_rec;
      output_rec.delimiter = " ";
      output_rec.setStream(ss);
      
      BOOST_CHECK_NO_THROW( output_rec << "x" << "2*x" << "x^2" );
      BOOST_CHECK_NO_THROW( output_rec << data_recorder::end_name_row );
      for(double x = 0; x < 10.1; x += 0.5) {
        BOOST_CHECK_NO_THROW( output_rec << x << 2*x << x*x );
        BOOST_CHECK_NO_THROW( output_rec << data_recorder::end_value_row );
      };
      BOOST_CHECK_NO_THROW( output_rec << data_recorder::flush );
    };
    
    {
      ascii_extractor input_rec;
      input_rec.delimiter = " ";
      input_rec.setStream(ss);
      
      BOOST_CHECK_EQUAL( input_rec.getColCount(), 3 );
      
      std::string s1, s2, s3;
      BOOST_CHECK_NO_THROW( input_rec >> s1 >> s2 >> s3 );
      BOOST_CHECK( s1 == "x" );
      BOOST_CHECK( s2 == "2*x" );
      BOOST_CHECK( s3 == "x^2" );
      for(double x = 0; x < 10.1; x += 0.5) {
        double v1, v2, v3;
        BOOST_CHECK_NO_THROW( input_rec >> v1 >> v2 >> v3 );
        BOOST_CHECK_CLOSE( v1, x, 1e-6 );
        BOOST_CHECK_CLOSE( v2, (2.0*x), 1e-6 );
        BOOST_CHECK_CLOSE( v3, (x*x), 1e-6 );
        BOOST_CHECK_NO_THROW( input_rec >> data_extractor::end_value_row );
      };
      BOOST_CHECK_NO_THROW( input_rec >> data_extractor::close );
    };
    
  };
  
};


BOOST_AUTO_TEST_CASE( ascii_tab_record_extract_test )
{
  using namespace ReaK;
  using namespace recorder;
  
  {
    std::stringstream ss;
    {
      ascii_recorder output_rec;
      output_rec.delimiter = "\t";
      output_rec.setStream(ss);
      
      BOOST_CHECK_NO_THROW( output_rec << "x" << "2*x" << "x^2" );
      BOOST_CHECK_NO_THROW( output_rec << data_recorder::end_name_row );
      for(double x = 0; x < 10.1; x += 0.5) {
        BOOST_CHECK_NO_THROW( output_rec << x << 2*x << x*x );
        BOOST_CHECK_NO_THROW( output_rec << data_recorder::end_value_row );
      };
      BOOST_CHECK_NO_THROW( output_rec << data_recorder::flush );
    };
    
    {
      ascii_extractor input_rec;
      input_rec.delimiter = "\t";
      input_rec.setStream(ss);
      
      BOOST_CHECK_EQUAL( input_rec.getColCount(), 3 );
      
      std::string s1, s2, s3;
      BOOST_CHECK_NO_THROW( input_rec >> s1 >> s2 >> s3 );
      BOOST_CHECK( s1 == "x" );
      BOOST_CHECK( s2 == "2*x" );
      BOOST_CHECK( s3 == "x^2" );
      for(double x = 0; x < 10.1; x += 0.5) {
        double v1, v2, v3;
        BOOST_CHECK_NO_THROW( input_rec >> v1 >> v2 >> v3 );
        BOOST_CHECK_CLOSE( v1, x, 1e-6 );
        BOOST_CHECK_CLOSE( v2, (2.0*x), 1e-6 );
        BOOST_CHECK_CLOSE( v3, (x*x), 1e-6 );
        BOOST_CHECK_NO_THROW( input_rec >> data_extractor::end_value_row );
      };
      BOOST_CHECK_NO_THROW( input_rec >> data_extractor::close );
    };
    
  };
  
};


BOOST_AUTO_TEST_CASE( ascii_comma_record_extract_test )
{
  using namespace ReaK;
  using namespace recorder;
  
  {
    std::stringstream ss;
    {
      ascii_recorder output_rec;
      output_rec.delimiter = ", ";
      output_rec.setStream(ss);
      
      BOOST_CHECK_NO_THROW( output_rec << "x" << "2*x" << "x^2" );
      BOOST_CHECK_NO_THROW( output_rec << data_recorder::end_name_row );
      for(double x = 0; x < 10.1; x += 0.5) {
        BOOST_CHECK_NO_THROW( output_rec << x << 2*x << x*x );
        BOOST_CHECK_NO_THROW( output_rec << data_recorder::end_value_row );
      };
      BOOST_CHECK_NO_THROW( output_rec << data_recorder::flush );
    };
    
    {
      ascii_extractor input_rec;
      input_rec.delimiter = ", ";
      input_rec.setStream(ss);
      
      BOOST_CHECK_EQUAL( input_rec.getColCount(), 3 );
      
      std::string s1, s2, s3;
      BOOST_CHECK_NO_THROW( input_rec >> s1 >> s2 >> s3 );
      BOOST_CHECK( s1 == "x" );
      BOOST_CHECK( s2 == "2*x" );
      BOOST_CHECK( s3 == "x^2" );
      for(double x = 0; x < 10.1; x += 0.5) {
        double v1, v2, v3;
        BOOST_CHECK_NO_THROW( input_rec >> v1 >> v2 >> v3 );
        BOOST_CHECK_CLOSE( v1, x, 1e-6 );
        BOOST_CHECK_CLOSE( v2, (2.0*x), 1e-6 );
        BOOST_CHECK_CLOSE( v3, (x*x), 1e-6 );
        BOOST_CHECK_NO_THROW( input_rec >> data_extractor::end_value_row );
      };
      BOOST_CHECK_NO_THROW( input_rec >> data_extractor::close );
    };
    
  };
  
};



// NOTE: also test the named-value-row input-output.
BOOST_AUTO_TEST_CASE( bin_record_extract_test )
{
  using namespace ReaK;
  using namespace recorder;
  
  {
    std::stringstream ss;
    {
      bin_recorder output_rec;
      output_rec.setStream(ss);
      
      BOOST_CHECK_NO_THROW( output_rec << "x" << "2*x" << "x^2" );
      BOOST_CHECK_NO_THROW( output_rec << data_recorder::end_name_row );
      named_value_row vr = output_rec.getFreshNamedValueRow();
      for(double x = 0; x < 10.1; x += 0.5) {
        vr["x"] = x;
        vr["2*x"] = 2*x;
        vr["x^2"] = x*x;
        BOOST_CHECK_NO_THROW( output_rec << vr );
      };
      BOOST_CHECK_NO_THROW( output_rec << data_recorder::flush );
    };
    
    {
      bin_extractor input_rec;
      input_rec.setStream(ss);
      
      BOOST_CHECK_EQUAL( input_rec.getColCount(), 3 );
      
      std::string s1, s2, s3;
      BOOST_CHECK_NO_THROW( input_rec >> s1 >> s2 >> s3 );
      BOOST_CHECK( s1 == "x" );
      BOOST_CHECK( s2 == "2*x" );
      BOOST_CHECK( s3 == "x^2" );
      named_value_row vr = input_rec.getFreshNamedValueRow();
      for(double x = 0; x < 10.1; x += 0.5) {
        BOOST_CHECK_NO_THROW( input_rec >> vr );
        BOOST_CHECK_CLOSE( vr["x"], x, 1e-6 );
        BOOST_CHECK_CLOSE( vr["2*x"], (2.0*x), 1e-6 );
        BOOST_CHECK_CLOSE( vr["x^2"], (x*x), 1e-6 );
      };
      BOOST_CHECK_NO_THROW( input_rec >> data_extractor::close );
    };
    
  };
  
};




struct net_server_runner {
  bool* succeeded;
  unsigned int* num_points;
  std::string server_uri;
  net_server_runner(bool* aSucceeded, unsigned int* aNumPoints, 
                    const std::string& aURI) : succeeded(aSucceeded), num_points(aNumPoints), server_uri(aURI) { };
  
  void operator()() {
    
    using namespace ReaK;
    using namespace recorder;
    
    try {
      network_recorder output_rec(server_uri);
      output_rec << "x" << "2*x" << "x^2" << data_recorder::end_name_row;
      for(double x = 0.0; x < 10.1; x += 0.5) {
        output_rec << x << 2*x << x*x << data_recorder::end_value_row;
        *num_points += 1;
      };
      output_rec << data_recorder::flush;
    } catch(...) {
      *succeeded = false;
      return;
    };
    *succeeded = true;
  };  
  
};


BOOST_AUTO_TEST_CASE( net_tcp_record_extract_test )
{
  
  using namespace ReaK;
  using namespace recorder;
  
  bool server_worked = false;
  unsigned int server_sent = 0;
  
  net_server_runner srv(&server_worked, &server_sent, "tcp:localhost:17020");
  ReaKaux::thread server_thd( srv );
  
  // it is necessary to give some time for the server to get up and waiting before a client can be created:
  ReaKaux::this_thread::sleep_for(ReaKaux::chrono::nanoseconds(10000));
  ReaKaux::this_thread::yield();
  
  network_extractor input_rec("tcp:localhost:17020");
  
  BOOST_CHECK_EQUAL( input_rec.getColCount(), 3 );
  std::string s1, s2, s3;
  BOOST_CHECK_NO_THROW( input_rec >> s1 >> s2 >> s3 );
  BOOST_CHECK( s1 == "x" );
  BOOST_CHECK( s2 == "2*x" );
  BOOST_CHECK( s3 == "x^2" );
  
  for(double x = 0; x < 10.1; x += 0.5) {
    double v1, v2, v3;
    BOOST_CHECK_NO_THROW( input_rec >> v1 >> v2 >> v3 );
    BOOST_CHECK_CLOSE( v1, x, 1e-6 );
    BOOST_CHECK_CLOSE( v2, (2.0*x), 1e-6 );
    BOOST_CHECK_CLOSE( v3, (x*x), 1e-6 );
    BOOST_CHECK_NO_THROW( input_rec >> data_extractor::end_value_row );
  };
  
  if( server_thd.joinable() )
    BOOST_CHECK_NO_THROW( server_thd.join() );
  BOOST_CHECK_EQUAL( server_sent, 21 );
  BOOST_CHECK( server_worked );
  
};


BOOST_AUTO_TEST_CASE( net_udp_record_extract_test )
{
  
  using namespace ReaK;
  using namespace recorder;
  
  bool server_worked = false;
  unsigned int server_sent = 0;
  
  net_server_runner srv(&server_worked, &server_sent, "udp:localhost:17021");
  ReaKaux::thread server_thd( srv );
  
  // it is necessary to give some time for the server to get up and waiting before a client can be created:
  ReaKaux::this_thread::sleep_for(ReaKaux::chrono::nanoseconds(10000));
  ReaKaux::this_thread::yield();
  
  network_extractor input_rec("udp:localhost:17021");
  
  BOOST_CHECK_EQUAL( input_rec.getColCount(), 3 );
  std::string s1, s2, s3;
  BOOST_CHECK_NO_THROW( input_rec >> s1 >> s2 >> s3 );
  BOOST_CHECK( s1 == "x" );
  BOOST_CHECK( s2 == "2*x" );
  BOOST_CHECK( s3 == "x^2" );
  
  for(double x = 0; x < 10.1; x += 0.5) {
    double v1, v2, v3;
    BOOST_CHECK_NO_THROW( input_rec >> v1 >> v2 >> v3 );
    BOOST_CHECK_CLOSE( v1, x, 1e-6 );
    BOOST_CHECK_CLOSE( v2, (2.0*x), 1e-6 );
    BOOST_CHECK_CLOSE( v3, (x*x), 1e-6 );
    BOOST_CHECK_NO_THROW( input_rec >> data_extractor::end_value_row );
  };
  
  if( server_thd.joinable() )
    BOOST_CHECK_NO_THROW( server_thd.join() );
  BOOST_CHECK_EQUAL( server_sent, 21 );
  BOOST_CHECK( server_worked );
  
};


BOOST_AUTO_TEST_CASE( net_raw_udp_record_extract_test )
{
  
  using namespace ReaK;
  using namespace recorder;
  
  bool server_worked = false;
  unsigned int server_sent = 0;
  
  network_extractor input_rec;
  input_rec.addName("x");
  input_rec.addName("2*x");
  input_rec.addName("x^2");
  
  input_rec.setFileName("raw_udp:localhost:17022");
  
  BOOST_CHECK_EQUAL( input_rec.getColCount(), 3 );
  std::string s1, s2, s3;
  BOOST_CHECK_NO_THROW( input_rec >> s1 >> s2 >> s3 );
  BOOST_CHECK( s1 == "x" );
  BOOST_CHECK( s2 == "2*x" );
  BOOST_CHECK( s3 == "x^2" );
  
  net_server_runner srv(&server_worked, &server_sent, "raw_udp:localhost:17022");
  ReaKaux::thread server_thd( srv );
  
  for(double x = 0; x < 10.1; x += 0.5) {
    double v1, v2, v3;
    BOOST_CHECK_NO_THROW( input_rec >> v1 >> v2 >> v3 );
    BOOST_CHECK_CLOSE( v1, x, 1e-6 );
    BOOST_CHECK_CLOSE( v2, (2.0*x), 1e-6 );
    BOOST_CHECK_CLOSE( v3, (x*x), 1e-6 );
    BOOST_CHECK_NO_THROW( input_rec >> data_extractor::end_value_row );
  };
  
  if( server_thd.joinable() )
    BOOST_CHECK_NO_THROW( server_thd.join() );
  BOOST_CHECK_EQUAL( server_sent, 21 );
  BOOST_CHECK( server_worked );
  
};







BOOST_AUTO_TEST_CASE( vector_record_extract_test )
{
  using namespace ReaK;
  using namespace recorder;
  
  {
    std::vector< std::vector<double> > vec;
    {
      vector_recorder output_rec;
      output_rec.setVecData(&vec);
      
      BOOST_CHECK_NO_THROW( output_rec << "x" << "2*x" << "x^2" );
      BOOST_CHECK_NO_THROW( output_rec << data_recorder::end_name_row );
      for(double x = 0; x < 10.1; x += 0.5) {
        BOOST_CHECK_NO_THROW( output_rec << x << 2*x << x*x );
        BOOST_CHECK_NO_THROW( output_rec << data_recorder::end_value_row );
      };
      BOOST_CHECK_NO_THROW( output_rec << data_recorder::flush );
    };
    
    {
      vector_extractor input_rec;
      input_rec.setVecData(&vec);
      
      input_rec.addName("x");
      input_rec.addName("2*x");
      input_rec.addName("x^2");
      BOOST_CHECK_EQUAL( input_rec.getColCount(), 3 );
      
      std::string s1, s2, s3;
      BOOST_CHECK_NO_THROW( input_rec >> s1 >> s2 >> s3 );
      BOOST_CHECK( s1 == "x" );
      BOOST_CHECK( s2 == "2*x" );
      BOOST_CHECK( s3 == "x^2" );
      for(double x = 0; x < 10.1; x += 0.5) {
        double v1, v2, v3;
        BOOST_CHECK_NO_THROW( input_rec >> v1 >> v2 >> v3 );
        BOOST_CHECK_CLOSE( v1, x, 1e-6 );
        BOOST_CHECK_CLOSE( v2, (2.0*x), 1e-6 );
        BOOST_CHECK_CLOSE( v3, (x*x), 1e-6 );
        BOOST_CHECK_NO_THROW( input_rec >> data_extractor::end_value_row );
      };
      BOOST_CHECK_NO_THROW( input_rec >> data_extractor::close );
    };
    
  };
  
};














