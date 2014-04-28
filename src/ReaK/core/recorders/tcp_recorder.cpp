
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

#include "tcp_recorder.hpp"

#include <boost/asio.hpp>



#include <stdint.h>

#ifdef WIN32

#include <winsock2.h>

#else

#include <netinet/in.h>
//#include <arpa/inet.h>

#endif

namespace ReaK {

namespace recorder {
  
namespace {

union double_to_ulong {
  double   d;
  uint64_t ui64;
  uint32_t ui32[2];
};


template <typename UnionT>
void ntoh_2ui32(UnionT& value) {
#if RK_BYTE_ORDER == RK_ORDER_LITTLE_ENDIAN
  uint32_t tmp = ntohl(value.ui32[0]);
  value.ui32[0] = ntohl(value.ui32[1]);
  value.ui32[1] = tmp;
#endif
  // NOTE: for 64-bit values, there is no point in supporting PDP-endianness, as 64-bit values are not supported by PDP platforms.
};

template <typename UnionT>
void hton_2ui32(UnionT& value) {
#if RK_BYTE_ORDER == RK_ORDER_LITTLE_ENDIAN
  uint32_t tmp = htonl(value.ui32[0]);
  value.ui32[0] = htonl(value.ui32[1]);
  value.ui32[1] = tmp;
#endif
  // NOTE: for 64-bit values, there is no point in supporting PDP-endianness, as 64-bit values are not supported by PDP platforms.
};


};



class tcp_server_impl {
  public:
    boost::asio::io_service io_service;
    boost::asio::ip::tcp::acceptor acceptor;
    boost::asio::ip::tcp::socket socket;
    boost::asio::basic_streambuf<> row_buf;
    
    tcp_server_impl(std::size_t port_num) : 
      io_service(), 
      acceptor(io_service, boost::asio::ip::tcp::endpoint(boost::asio::ip::tcp::v4(), port_num)),
      socket(io_service) {
      acceptor.accept(socket);
    };
    
};

class tcp_client_impl {
  public:
    boost::asio::io_service io_service;
    boost::asio::ip::tcp::endpoint endpoint;
    boost::asio::ip::tcp::socket socket;
    boost::asio::basic_streambuf<> row_buf;
    
    tcp_client_impl(const std::string& ip4_address, std::size_t port_num) : 
      io_service(), 
      endpoint(boost::asio::ip::address_v4::from_string(ip4_address), port_num),
      socket(io_service) { 
      socket.connect(endpoint);
    };
    
};



tcp_recorder::tcp_recorder() : data_recorder(), pimpl(), apply_network_order(false) { };

tcp_recorder::tcp_recorder(const std::string& aFileName) {
  setFileName(aFileName);
};

tcp_recorder::~tcp_recorder() { };

void tcp_recorder::writeRow() {
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  if((pimpl) && (pimpl->socket.is_open()) && (rowCount > 0) && (colCount > 0)) {
    std::ostream s_tmp(&(pimpl->row_buf));
    {
      for(std::size_t i = 0; i < colCount; ++i) {
        if(apply_network_order) {
          double_to_ulong tmp; tmp.d = values_rm.front();
          hton_2ui32(tmp);
          s_tmp.write(reinterpret_cast<char*>(&tmp),sizeof(double));
        } else {
          double tmp(values_rm.front());
          s_tmp.write(reinterpret_cast<char*>(&tmp),sizeof(double));
        };
        values_rm.pop();
      };
      --rowCount;
    };
    std::size_t len = boost::asio::write(pimpl->socket, pimpl->row_buf);
    pimpl->row_buf.consume(len);
  };
};

void tcp_recorder::writeNames() {
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  if((pimpl) && (pimpl->socket.is_open())) {
    std::stringstream ss;
    for(std::vector<std::string>::iterator it = names.begin(); it != names.end(); ++it)
      ss << " " << *it;
    std::string data_str = ss.str();
    uint32_t data_len = htonl(data_str.size());
    std::ostream s_tmp(&(pimpl->row_buf));
    s_tmp.write(reinterpret_cast<char*>(&data_len), sizeof(uint32_t));
    s_tmp.write(data_str.c_str(), data_str.size());
    std::size_t len = boost::asio::write(pimpl->socket, pimpl->row_buf);
    pimpl->row_buf.consume(len);
  };
};

void tcp_recorder::setFileName(const std::string& aFileName) {
  if(colCount != 0) {
    *this << close;
    
    ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
    std::size_t port_num = 0;
    std::stringstream ss(aFileName); ss >> port_num;
    pimpl = shared_ptr<tcp_server_impl>(new tcp_server_impl(port_num));
    colCount = names.size();
    lock_here.unlock();
    writeNames();
  } else {
    ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
    std::size_t port_num = 0;
    std::stringstream ss(aFileName); ss >> port_num;
    pimpl = shared_ptr<tcp_server_impl>(new tcp_server_impl(port_num));
  };
};





tcp_extractor::tcp_extractor() : data_extractor(), pimpl(), apply_network_order(false) { };

tcp_extractor::tcp_extractor(const std::string& aFileName) : data_extractor(), pimpl() {
  setFileName(aFileName);
};

tcp_extractor::~tcp_extractor() {};

bool tcp_extractor::readRow() {
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  shared_ptr<tcp_client_impl> pimpl_tmp = pimpl;
  if((pimpl_tmp) && (pimpl_tmp->socket.is_open()) && (colCount > 0)) {
    try {
      boost::asio::streambuf::mutable_buffers_type bufs = pimpl_tmp->row_buf.prepare(colCount * sizeof(double));
      std::size_t len = boost::asio::read(pimpl_tmp->socket, bufs);
      pimpl_tmp->row_buf.commit(len);
      if(len < colCount * sizeof(double))
        return false;
    } catch(...) {
      return false;
    };
    std::istream s_tmp(&(pimpl_tmp->row_buf));
    for(std::size_t i = 0; (i < colCount) && (s_tmp); ++i) {
      if(apply_network_order) {
        double_to_ulong tmp;
        s_tmp.read(reinterpret_cast<char*>(&tmp),sizeof(double));
        ntoh_2ui32(tmp);
        values_rm.push(tmp.d);
      } else {
        double tmp = 0;
        s_tmp.read(reinterpret_cast<char*>(&tmp),sizeof(double));
        values_rm.push(tmp);
      };
    };
  };
  return true;
};

bool tcp_extractor::readNames() {
  shared_ptr<tcp_client_impl> pimpl_tmp = pimpl;
  if((pimpl_tmp) && (pimpl_tmp->socket.is_open())) {
    uint32_t data_len = 0;
    {
      boost::asio::streambuf::mutable_buffers_type bufs = pimpl_tmp->row_buf.prepare(sizeof(uint32_t));
      std::size_t len = boost::asio::read(pimpl_tmp->socket, bufs);
      pimpl_tmp->row_buf.commit(len);
      std::istream s_tmp(&(pimpl_tmp->row_buf));
      s_tmp.read(reinterpret_cast<char*>(&data_len),sizeof(uint32_t));
      data_len = ntohl(data_len);
    };
    boost::asio::streambuf::mutable_buffers_type bufs = pimpl_tmp->row_buf.prepare(data_len);
    std::size_t len = boost::asio::read(pimpl_tmp->socket, bufs);
    pimpl_tmp->row_buf.commit(len);
    std::istream s_tmp(&(pimpl_tmp->row_buf));
    std::string tmp_name = "";
    while(s_tmp >> tmp_name) {
      names.push_back(tmp_name);
      ++colCount;
    };
  };
  return true;
};

void tcp_extractor::setFileName(const std::string& aFileName) {
  if(colCount != 0)
    *this >> close;

  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  std::size_t i = aFileName.find(':');
  std::string ip4addr = aFileName.substr(0, i);
  std::size_t portnum = 17000;
  if(++i < aFileName.size()) {
    std::stringstream ss( aFileName.substr(i, aFileName.size() - i) );
    ss >> portnum;
  };
  pimpl = shared_ptr<tcp_client_impl>(new tcp_client_impl(ip4addr, portnum));
  readNames();
};



};


};



