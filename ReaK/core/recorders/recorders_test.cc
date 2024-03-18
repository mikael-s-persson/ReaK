
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

#include "ReaK/core/recorders/ascii_recorder.h"
#include "ReaK/core/recorders/bin_recorder.h"
#include "ReaK/core/recorders/vector_recorder.h"

#ifdef ENABLE_NETWORK_RECORDER
#include "ReaK/core/recorders/network_recorder.h"
#endif // ENABLE_NETWORK_RECORDER

#include <sstream>

#include <chrono>
#include <thread>

#include "gtest/gtest.h"

namespace ReaK::recorder {
namespace {

TEST(RecordersTests, AsciiSpaceRecordExtract) {
  {
    std::stringstream ss;
    {
      ascii_recorder output_rec;
      output_rec.delimiter = " ";
      output_rec.setStream(ss);

      EXPECT_NO_THROW(output_rec << "x"
                                 << "2*x"
                                 << "x^2");
      EXPECT_NO_THROW(output_rec << data_recorder::end_name_row);
      for (double x = 0; x < 10.1; x += 0.5) {
        EXPECT_NO_THROW(output_rec << x << 2 * x << x * x);
        EXPECT_NO_THROW(output_rec << data_recorder::end_value_row);
      }
      EXPECT_NO_THROW(output_rec << data_recorder::flush);
    }

    {
      ascii_extractor input_rec;
      input_rec.delimiter = " ";
      input_rec.setStream(ss);

      EXPECT_EQ(input_rec.getColCount(), 3);

      std::string s1;
      std::string s2;
      std::string s3;
      EXPECT_NO_THROW(input_rec >> s1 >> s2 >> s3);
      EXPECT_EQ(s1, "x");
      EXPECT_EQ(s2, "2*x");
      EXPECT_EQ(s3, "x^2");
      for (double x = 0; x < 10.1; x += 0.5) {
        double v1 = 0.0;
        double v2 = 0.0;
        double v3 = 0.0;
        EXPECT_NO_THROW(input_rec >> v1 >> v2 >> v3);
        EXPECT_NEAR(v1, x, 1e-6);
        EXPECT_NEAR(v2, (2.0 * x), 1e-6);
        EXPECT_NEAR(v3, (x * x), 1e-6);
        EXPECT_NO_THROW(input_rec >> data_extractor::end_value_row);
      }
      EXPECT_NO_THROW(input_rec >> data_extractor::close);
    }
  }
}

TEST(RecordersTests, AsciiTabRecordExtract) {
  {
    std::stringstream ss;
    {
      ascii_recorder output_rec;
      output_rec.delimiter = "\t";
      output_rec.setStream(ss);

      EXPECT_NO_THROW(output_rec << "x"
                                 << "2*x"
                                 << "x^2");
      EXPECT_NO_THROW(output_rec << data_recorder::end_name_row);
      for (double x = 0; x < 10.1; x += 0.5) {
        EXPECT_NO_THROW(output_rec << x << 2 * x << x * x);
        EXPECT_NO_THROW(output_rec << data_recorder::end_value_row);
      }
      EXPECT_NO_THROW(output_rec << data_recorder::flush);
    }

    {
      ascii_extractor input_rec;
      input_rec.delimiter = "\t";
      input_rec.setStream(ss);

      EXPECT_EQ(input_rec.getColCount(), 3);

      std::string s1;
      std::string s2;
      std::string s3;
      EXPECT_NO_THROW(input_rec >> s1 >> s2 >> s3);
      EXPECT_EQ(s1, "x");
      EXPECT_EQ(s2, "2*x");
      EXPECT_EQ(s3, "x^2");
      for (double x = 0; x < 10.1; x += 0.5) {
        double v1 = 0.0;
        double v2 = 0.0;
        double v3 = 0.0;
        EXPECT_NO_THROW(input_rec >> v1 >> v2 >> v3);
        EXPECT_NEAR(v1, x, 1e-6);
        EXPECT_NEAR(v2, (2.0 * x), 1e-6);
        EXPECT_NEAR(v3, (x * x), 1e-6);
        EXPECT_NO_THROW(input_rec >> data_extractor::end_value_row);
      }
      EXPECT_NO_THROW(input_rec >> data_extractor::close);
    }
  }
}

TEST(RecordersTests, AsciiCommaRecordExtract) {
  {
    std::stringstream ss;
    {
      ascii_recorder output_rec;
      output_rec.delimiter = ", ";
      output_rec.setStream(ss);

      EXPECT_NO_THROW(output_rec << "x"
                                 << "2*x"
                                 << "x^2");
      EXPECT_NO_THROW(output_rec << data_recorder::end_name_row);
      for (double x = 0; x < 10.1; x += 0.5) {
        EXPECT_NO_THROW(output_rec << x << 2 * x << x * x);
        EXPECT_NO_THROW(output_rec << data_recorder::end_value_row);
      }
      EXPECT_NO_THROW(output_rec << data_recorder::flush);
    }

    {
      ascii_extractor input_rec;
      input_rec.delimiter = ", ";
      input_rec.setStream(ss);

      EXPECT_EQ(input_rec.getColCount(), 3);

      std::string s1;
      std::string s2;
      std::string s3;
      EXPECT_NO_THROW(input_rec >> s1 >> s2 >> s3);
      EXPECT_EQ(s1, "x");
      EXPECT_EQ(s2, "2*x");
      EXPECT_EQ(s3, "x^2");
      for (double x = 0; x < 10.1; x += 0.5) {
        double v1 = 0.0;
        double v2 = 0.0;
        double v3 = 0.0;
        EXPECT_NO_THROW(input_rec >> v1 >> v2 >> v3);
        EXPECT_NEAR(v1, x, 1e-6);
        EXPECT_NEAR(v2, (2.0 * x), 1e-6);
        EXPECT_NEAR(v3, (x * x), 1e-6);
        EXPECT_NO_THROW(input_rec >> data_extractor::end_value_row);
      }
      EXPECT_NO_THROW(input_rec >> data_extractor::close);
    }
  }
}

// NOTE: also test the named-value-row input-output.
TEST(RecordersTests, BinRecordExtract) {
  {
    std::stringstream ss;
    {
      bin_recorder output_rec;
      output_rec.setStream(ss);

      EXPECT_NO_THROW(output_rec << "x"
                                 << "2*x"
                                 << "x^2");
      EXPECT_NO_THROW(output_rec << data_recorder::end_name_row);
      named_value_row vr = output_rec.getFreshNamedValueRow();
      for (double x = 0; x < 10.1; x += 0.5) {
        vr["x"] = x;
        vr["2*x"] = 2 * x;
        vr["x^2"] = x * x;
        EXPECT_NO_THROW(output_rec << vr);
      }
      EXPECT_NO_THROW(output_rec << data_recorder::flush);
    }

    {
      bin_extractor input_rec;
      input_rec.setStream(ss);

      EXPECT_EQ(input_rec.getColCount(), 3);

      std::string s1;
      std::string s2;
      std::string s3;
      EXPECT_NO_THROW(input_rec >> s1 >> s2 >> s3);
      EXPECT_EQ(s1, "x");
      EXPECT_EQ(s2, "2*x");
      EXPECT_EQ(s3, "x^2");
      named_value_row vr = input_rec.getFreshNamedValueRow();
      for (double x = 0; x < 10.1; x += 0.5) {
        EXPECT_NO_THROW(input_rec >> vr);
        EXPECT_NEAR(vr["x"], x, 1e-6);
        EXPECT_NEAR(vr["2*x"], (2.0 * x), 1e-6);
        EXPECT_NEAR(vr["x^2"], (x * x), 1e-6);
      }
      EXPECT_NO_THROW(input_rec >> data_extractor::close);
    }
  }
}

#ifdef ENABLE_NETWORK_RECORDER
struct net_server_runner {
  bool* succeeded;
  unsigned int* num_points;
  std::string server_uri;
  net_server_runner(bool* aSucceeded, unsigned int* aNumPoints,
                    const std::string& aURI)
      : succeeded(aSucceeded), num_points(aNumPoints), server_uri(aURI) {}

  void operator()() const {

    using namespace ReaK;
    using namespace recorder;

    try {
      network_recorder output_rec(server_uri);
      output_rec << "x"
                 << "2*x"
                 << "x^2" << data_recorder::end_name_row;
      for (double x = 0.0; x < 10.1; x += 0.5) {
        output_rec << x << 2 * x << x * x << data_recorder::end_value_row;
        *num_points += 1;
      }
      output_rec << data_recorder::flush;
    } catch (...) {
      *succeeded = false;
      return;
    }
    *succeeded = true;
  }
};

TEST(RecordersTests, NetTcpRecordExtract) {
  bool server_worked = false;
  unsigned int server_sent = 0;

  net_server_runner srv(&server_worked, &server_sent, "tcp:localhost:17020");
  std::thread server_thd(srv);

  // it is necessary to give some time for the server to get up and waiting before a client can be created:
  std::this_thread::sleep_for(std::chrono::nanoseconds(10000));
  std::this_thread::yield();

  network_extractor input_rec("tcp:localhost:17020");

  EXPECT_EQ(input_rec.getColCount(), 3);
  std::string s1;
  std::string s2;
  std::string s3;
  EXPECT_NO_THROW(input_rec >> s1 >> s2 >> s3);
  EXPECT_EQ(s1, "x");
  EXPECT_EQ(s2, "2*x");
  EXPECT_EQ(s3, "x^2");

  for (double x = 0; x < 10.1; x += 0.5) {
    double v1 = 0.0;
    double v2 = 0.0;
    double v3 = 0.0;
    EXPECT_NO_THROW(input_rec >> v1 >> v2 >> v3);
    EXPECT_NEAR(v1, x, 1e-6);
    EXPECT_NEAR(v2, (2.0 * x), 1e-6);
    EXPECT_NEAR(v3, (x * x), 1e-6);
    EXPECT_NO_THROW(input_rec >> data_extractor::end_value_row);
  }

  if (server_thd.joinable()) {
    EXPECT_NO_THROW(server_thd.join());
  }
  EXPECT_EQ(server_sent, 21);
  EXPECT_TRUE(server_worked);
}

TEST(RecordersTests, NetUdpRecordExtract) {
  bool server_worked = false;
  unsigned int server_sent = 0;

  net_server_runner srv(&server_worked, &server_sent, "udp:localhost:17021");
  std::thread server_thd(srv);

  // it is necessary to give some time for the server to get up and waiting before a client can be created:
  std::this_thread::sleep_for(std::chrono::nanoseconds(10000));
  std::this_thread::yield();

  network_extractor input_rec("udp:localhost:17021");

  EXPECT_EQ(input_rec.getColCount(), 3);
  std::string s1;
  std::string s2;
  std::string s3;
  EXPECT_NO_THROW(input_rec >> s1 >> s2 >> s3);
  EXPECT_EQ(s1, "x");
  EXPECT_EQ(s2, "2*x");
  EXPECT_EQ(s3, "x^2");

  for (double x = 0; x < 10.1; x += 0.5) {
    double v1 = 0.0;
    double v2 = 0.0;
    double v3 = 0.0;
    EXPECT_NO_THROW(input_rec >> v1 >> v2 >> v3);
    EXPECT_NEAR(v1, x, 1e-6);
    EXPECT_NEAR(v2, (2.0 * x), 1e-6);
    EXPECT_NEAR(v3, (x * x), 1e-6);
    EXPECT_NO_THROW(input_rec >> data_extractor::end_value_row);
  }

  if (server_thd.joinable()) {
    EXPECT_NO_THROW(server_thd.join());
  }
  EXPECT_EQ(server_sent, 21);
  EXPECT_TRUE(server_worked);
}

TEST(RecordersTests, NetRawUdpRecordExtract) {
  bool server_worked = false;
  unsigned int server_sent = 0;

  network_extractor input_rec;
  input_rec.addName("x");
  input_rec.addName("2*x");
  input_rec.addName("x^2");

  input_rec.setFileName("raw_udp:localhost:17022");

  EXPECT_EQ(input_rec.getColCount(), 3);
  std::string s1;
  std::string s2;
  std::string s3;
  EXPECT_NO_THROW(input_rec >> s1 >> s2 >> s3);
  EXPECT_EQ(s1, "x");
  EXPECT_EQ(s2, "2*x");
  EXPECT_EQ(s3, "x^2");

  net_server_runner srv(&server_worked, &server_sent,
                        "raw_udp:localhost:17022");
  std::thread server_thd(srv);

  for (double x = 0; x < 10.1; x += 0.5) {
    double v1 = 0.0;
    double v2 = 0.0;
    double v3 = 0.0;
    EXPECT_NO_THROW(input_rec >> v1 >> v2 >> v3);
    EXPECT_NEAR(v1, x, 1e-6);
    EXPECT_NEAR(v2, (2.0 * x), 1e-6);
    EXPECT_NEAR(v3, (x * x), 1e-6);
    EXPECT_NO_THROW(input_rec >> data_extractor::end_value_row);
  }

  if (server_thd.joinable()) {
    EXPECT_NO_THROW(server_thd.join());
  }
  EXPECT_EQ(server_sent, 21);
  EXPECT_TRUE(server_worked);
}
#endif // ENABLE_NETWORK_RECORDER

TEST(RecordersTests, VectorRecordExtract) {
  {
    std::vector<std::vector<double>> vec;
    {
      vector_recorder output_rec;
      output_rec.setVecData(&vec);

      EXPECT_NO_THROW(output_rec << "x"
                                 << "2*x"
                                 << "x^2");
      EXPECT_NO_THROW(output_rec << data_recorder::end_name_row);
      for (double x = 0; x < 10.1; x += 0.5) {
        EXPECT_NO_THROW(output_rec << x << 2 * x << x * x);
        EXPECT_NO_THROW(output_rec << data_recorder::end_value_row);
      }
      EXPECT_NO_THROW(output_rec << data_recorder::flush);
    }

    {
      vector_extractor input_rec;
      input_rec.setVecData(&vec);

      input_rec.addName("x");
      input_rec.addName("2*x");
      input_rec.addName("x^2");
      EXPECT_EQ(input_rec.getColCount(), 3);

      std::string s1;
      std::string s2;
      std::string s3;
      EXPECT_NO_THROW(input_rec >> s1 >> s2 >> s3);
      EXPECT_EQ(s1, "x");
      EXPECT_EQ(s2, "2*x");
      EXPECT_EQ(s3, "x^2");
      for (double x = 0; x < 10.1; x += 0.5) {
        double v1 = 0.0;
        double v2 = 0.0;
        double v3 = 0.0;
        EXPECT_NO_THROW(input_rec >> v1 >> v2 >> v3);
        EXPECT_NEAR(v1, x, 1e-6);
        EXPECT_NEAR(v2, (2.0 * x), 1e-6);
        EXPECT_NEAR(v3, (x * x), 1e-6);
        EXPECT_NO_THROW(input_rec >> data_extractor::end_value_row);
      }
      EXPECT_NO_THROW(input_rec >> data_extractor::close);
    }
  }
}

}  // namespace
}  // namespace ReaK::recorder
