
/*
 *    Copyright 2014 Sven Mikael Persson
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

#include <ReaK/core/recorders/network_recorder.hpp>

#include <ReaK/core/base/endian_conversions.hpp>

#include <boost/asio.hpp>


namespace ReaK {

namespace recorder {



class network_server_impl {
  public:
    boost::asio::basic_streambuf<> row_buf;
    
    network_server_impl() { };
    virtual ~network_server_impl() { };
    
    virtual bool isOpen() const = 0;
    virtual std::size_t writeRowBuffer() = 0;
    
    virtual void writeRow(std::queue<double>& values_rm, unsigned int colCount) {
      std::ostream s_tmp(&row_buf);
      for(std::size_t i = 0; i < colCount; ++i) {
        double_to_ulong tmp; tmp.d = values_rm.front();
        hton_2ui32(tmp);
        s_tmp.write(reinterpret_cast<char*>(&tmp),sizeof(double));
        values_rm.pop();
      };
      writeRowBuffer();
    };
    
    virtual void writeNames(const std::vector<std::string>& names) {
      std::stringstream ss;
      for(std::vector<std::string>::const_iterator it = names.begin(); it != names.end(); ++it)
        ss << " " << *it;
      std::string data_str = ss.str();
      uint32_t data_len = htonl(data_str.size());
      std::ostream s_tmp(&row_buf);
      s_tmp.write(reinterpret_cast<char*>(&data_len), sizeof(uint32_t));
      writeRowBuffer();  // <-- this was done in original UDP server.
      s_tmp.write(data_str.c_str(), data_str.size());
      writeRowBuffer();
    };
    
};

class network_client_impl {
  public:
    boost::asio::basic_streambuf<> row_buf;
    
    network_client_impl() { };
    virtual ~network_client_impl() { };
    
    virtual bool isOpen() const = 0;
    virtual std::size_t readRowBuffer(std::size_t data_len) = 0;
    
    virtual bool readRow(std::queue<double>& values_rm, unsigned int colCount) {
      try {
        std::size_t len = readRowBuffer(colCount * sizeof(double));
        if(len == std::numeric_limits< std::size_t >::max())
          return true;
        if(len < colCount * sizeof(double))
          return false;
      } catch(...) {
        return false;
      };
      std::istream s_tmp(&row_buf);
      for(std::size_t i = 0; (i < colCount) && (s_tmp); ++i) {
        double_to_ulong tmp;
        s_tmp.read(reinterpret_cast<char*>(&tmp),sizeof(double));
        ntoh_2ui32(tmp);
        values_rm.push(tmp.d);
//         RK_NOTICE(1," Received: " << values_rm.back());
      };
      return true;
    };
    
    virtual void readNames(std::vector<std::string>& names) {
      names.clear();
      uint32_t data_len = 0;
      {
        while(readRowBuffer(sizeof(uint32_t)) == std::numeric_limits< std::size_t >::max())
          ;
        std::istream s_tmp(&row_buf);
        s_tmp.read(reinterpret_cast<char*>(&data_len),sizeof(uint32_t));
        data_len = ntohl(data_len);
      };
      while(readRowBuffer(data_len) == std::numeric_limits< std::size_t >::max())
        ;
      std::istream s_tmp(&row_buf);
      std::string tmp_name = "";
      while(s_tmp >> tmp_name)
        names.push_back(tmp_name);
    };
};


namespace detail {

class tcp_server_impl : public network_server_impl {
  public:
    boost::asio::io_service io_service;
    boost::asio::ip::tcp::acceptor acceptor;
    boost::asio::ip::tcp::socket socket;
    
    tcp_server_impl(std::size_t port_num) : 
      io_service(), 
      acceptor(io_service, boost::asio::ip::tcp::endpoint(boost::asio::ip::tcp::v4(), port_num)),
      socket(io_service) {
      acceptor.accept(socket);
    };
    
    bool isOpen() const { return socket.is_open(); };
    
    std::size_t writeRowBuffer() {
      std::size_t len = boost::asio::write(socket, row_buf);
      row_buf.consume(len);
      return len;
    };
    
};


class tcp_client_impl : public network_client_impl {
  public:
    boost::asio::io_service io_service;
    boost::asio::ip::tcp::endpoint endpoint;
    boost::asio::ip::tcp::socket socket;
    
    tcp_client_impl(const std::string& ip4_address, std::size_t port_num) : 
      io_service(), 
      endpoint(),
      socket(io_service) {
      
      boost::asio::ip::tcp::resolver addr_resolver(io_service);
      boost::asio::ip::tcp::resolver::query addr_query(boost::asio::ip::tcp::v4(), ip4_address, "");
      boost::asio::ip::tcp::resolver::iterator it = addr_resolver.resolve(addr_query);
      if(it != boost::asio::ip::tcp::resolver::iterator())
        endpoint = boost::asio::ip::tcp::endpoint(it->endpoint().address(), port_num);
      else
        throw std::invalid_argument("Could not resolve the host-name specified to a valid IPv4 address!");
      
      socket.connect(endpoint);
    };
    
    bool isOpen() const { return socket.is_open(); };
    
    std::size_t readRowBuffer(std::size_t data_len) {
      boost::asio::streambuf::mutable_buffers_type bufs = row_buf.prepare(data_len);
      std::size_t len = boost::asio::read(socket, bufs);
      row_buf.commit(len);
      return len;
    };
    
};



class udp_server_impl : public network_server_impl {
  public:
    boost::asio::io_service io_service;
    boost::asio::ip::udp::endpoint endpoint;
    boost::asio::ip::udp::socket socket;
    
    udp_server_impl(std::size_t port_num) : 
      io_service(), 
      endpoint(boost::asio::ip::udp::v4(), port_num), 
      socket(io_service) {
      
      {
        boost::asio::io_service tcp_io_service;
        boost::asio::ip::tcp::acceptor tcp_acceptor(tcp_io_service, boost::asio::ip::tcp::endpoint(boost::asio::ip::tcp::v4(), port_num));
        boost::asio::ip::tcp::socket tcp_socket(tcp_io_service);
        tcp_acceptor.accept(tcp_socket);
        endpoint.address(tcp_socket.remote_endpoint().address());
      };
      
      socket.open(boost::asio::ip::udp::v4());
      socket.set_option(boost::asio::ip::udp::socket::reuse_address(true));
    };
    
    bool isOpen() const { return socket.is_open(); };
    
    std::size_t writeRowBuffer() {
      std::size_t len = socket.send_to(row_buf.data(), endpoint);
      row_buf.consume(len);
      return len;
    };
    
};

class udp_client_impl : public network_client_impl {
  public:
    boost::asio::io_service io_service;
    boost::asio::ip::udp::endpoint endpoint;
    boost::asio::ip::udp::socket socket;
    
    udp_client_impl(const std::string& ip4_address, std::size_t port_num) : 
      io_service(), 
      endpoint(), 
      socket(io_service) { 
      
      boost::asio::ip::udp::resolver addr_resolver(io_service);
      boost::asio::ip::udp::resolver::query addr_query(boost::asio::ip::udp::v4(), ip4_address, "");
      boost::asio::ip::udp::resolver::iterator it = addr_resolver.resolve(addr_query);
      if(it != boost::asio::ip::udp::resolver::iterator())
        endpoint = boost::asio::ip::udp::endpoint(it->endpoint().address(), port_num);
      else
        throw std::invalid_argument("Could not resolve the host-name specified to a valid IPv4 address!");
      
      socket.open(boost::asio::ip::udp::v4());
      socket.set_option(boost::asio::ip::udp::socket::reuse_address(true));
      socket.bind(endpoint);
      
      {
        boost::asio::io_service tcp_io_service;
        boost::asio::ip::tcp::endpoint tcp_endpoint(endpoint.address(), endpoint.port());
        boost::asio::ip::tcp::socket tcp_socket(tcp_io_service);
        tcp_socket.connect(tcp_endpoint);
      };
      
    };
    
    bool isOpen() const { return socket.is_open(); };
    
    std::size_t readRowBuffer(std::size_t data_len) {
      if(socket.available() < data_len)
        return std::numeric_limits< std::size_t >::max();
      boost::asio::streambuf::mutable_buffers_type bufs = row_buf.prepare(data_len);
      std::size_t len = socket.receive_from(bufs, endpoint);
      row_buf.commit(len);
      return len;
    };
    
};


class raw_udp_server_impl : public network_server_impl {
  public:
    boost::asio::io_service io_service;
    boost::asio::ip::udp::endpoint endpoint;
    boost::asio::ip::udp::socket socket;
    
    raw_udp_server_impl(const std::string& ip4_address, std::size_t port_num) : 
      io_service(), 
      endpoint(), 
      socket(io_service) {
      
      boost::asio::ip::udp::resolver addr_resolver(io_service);
      boost::asio::ip::udp::resolver::query addr_query(boost::asio::ip::udp::v4(), ip4_address, "");
      boost::asio::ip::udp::resolver::iterator it = addr_resolver.resolve(addr_query);
      if(it != boost::asio::ip::udp::resolver::iterator())
        endpoint = boost::asio::ip::udp::endpoint(it->endpoint().address(), port_num);
      else
        throw std::invalid_argument("Could not resolve the host-name specified to a valid IPv4 address!");
      
      socket.open(boost::asio::ip::udp::v4());
      socket.set_option(boost::asio::ip::udp::socket::reuse_address(true));
    };
    
    bool isOpen() const { return socket.is_open(); };
    
    std::size_t writeRowBuffer() {
      std::size_t len = socket.send_to(row_buf.data(), endpoint);
      row_buf.consume(len);
      return len;
    };
    
    void writeNames(const std::vector<std::string>& names) { };
};

class raw_udp_client_impl : public network_client_impl {
  public:
    boost::asio::io_service io_service;
    boost::asio::ip::udp::endpoint endpoint;
    boost::asio::ip::udp::socket socket;
    
    raw_udp_client_impl(const std::string& ip4_address, std::size_t port_num) : 
      io_service(), 
      endpoint(), 
      socket(io_service) {
      
      boost::asio::ip::udp::resolver addr_resolver(io_service);
      boost::asio::ip::udp::resolver::query addr_query(boost::asio::ip::udp::v4(), ip4_address, "");
      boost::asio::ip::udp::resolver::iterator it = addr_resolver.resolve(addr_query);
      if(it != boost::asio::ip::udp::resolver::iterator())
        endpoint = boost::asio::ip::udp::endpoint(it->endpoint().address(), port_num);
      else
        throw std::invalid_argument("Could not resolve the host-name specified to a valid IPv4 address!");
      
      socket.open(boost::asio::ip::udp::v4());
      socket.set_option(boost::asio::ip::udp::socket::reuse_address(true));
      socket.bind(endpoint);
    };
    
    bool isOpen() const { return socket.is_open(); };
    
    std::size_t readRowBuffer(std::size_t data_len) {
      if(socket.available() < data_len)
        return std::numeric_limits< std::size_t >::max();
      boost::asio::streambuf::mutable_buffers_type bufs = row_buf.prepare(data_len);
      std::size_t len = socket.receive_from(bufs, endpoint);
      row_buf.commit(len);
      return len;
    };
    
    void readNames(std::vector<std::string>&) { };
};

};










network_recorder::network_recorder() : data_recorder(), pimpl() { };

network_recorder::network_recorder(const std::string& aFileName) {
  setFileName(aFileName);
};

network_recorder::~network_recorder() {
  closeRecordProcess();
};

void network_recorder::writeRow() {
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  shared_ptr<network_server_impl> pimpl_tmp = pimpl;
  if((pimpl_tmp) && (pimpl_tmp->isOpen()) && (rowCount > 0) && (colCount > 0)) {
    pimpl_tmp->writeRow(values_rm, names.size());
    --rowCount;
  };
};

void network_recorder::writeNames() {
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  shared_ptr<network_server_impl> pimpl_tmp = pimpl;
  if((pimpl_tmp) && (pimpl_tmp->isOpen()) && (colCount > 0))
    pimpl_tmp->writeNames(names);
};

void network_recorder::setFileName(const std::string& aFileName) {
  if(colCount != 0)
    *this << close;
    
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  
  std::string proto   = "tcp";
  std::string ip4addr = "localhost";
  std::size_t portnum = 17000;
  std::stringstream ss(aFileName);
  std::getline(ss, proto, ':');
  std::getline(ss, ip4addr, ':');
  ss >> portnum;
  
  if(proto == "tcp")
    pimpl = shared_ptr<detail::tcp_server_impl>(new detail::tcp_server_impl(portnum));
  else if(proto == "udp")
    pimpl = shared_ptr<detail::udp_server_impl>(new detail::udp_server_impl(portnum));
  else if(proto == "raw_udp")
    pimpl = shared_ptr<detail::raw_udp_server_impl>(new detail::raw_udp_server_impl(ip4addr, portnum));
  
  colCount = names.size();
  lock_here.unlock();
  writeNames();
};




network_extractor::network_extractor() : data_extractor(), pimpl() { };

network_extractor::network_extractor(const std::string& aFileName) : data_extractor(), pimpl() {
  setFileName(aFileName);
};

network_extractor::~network_extractor() {
  closeExtractProcess();
};

void network_extractor::addName(const std::string& s) { 
  if(colCount != 0)
    *this >> close;
  
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  names.push_back(s); 
};

bool network_extractor::readRow() {
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  shared_ptr<network_client_impl> pimpl_tmp = pimpl;
  if((pimpl_tmp) && (pimpl_tmp->isOpen()) && (colCount > 0))
    return pimpl_tmp->readRow(values_rm, names.size());
  return true;
};

bool network_extractor::readNames() {
  shared_ptr<network_client_impl> pimpl_tmp = pimpl;
  if((pimpl_tmp) && (pimpl_tmp->isOpen())) {
    pimpl_tmp->readNames(names);
    colCount = names.size();
  };
  return true;
};

void network_extractor::setFileName(const std::string& aFileName) {
  if(colCount != 0)
    *this >> close;
  
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  
  // Example filename (or URI): "tcp:localhost:17000" or "raw_udp:192.168.0.42:16069"
  std::string proto   = "tcp";
  std::string ip4addr = "localhost";
  std::size_t portnum = 17000;
  std::stringstream ss(aFileName);
  std::getline(ss, proto, ':');
  std::getline(ss, ip4addr, ':');
  ss >> portnum;
  
  if(proto == "tcp") {
    pimpl = shared_ptr<detail::tcp_client_impl>(new detail::tcp_client_impl(ip4addr, portnum));
  } else if(proto == "udp") {
    pimpl = shared_ptr<detail::udp_client_impl>(new detail::udp_client_impl(ip4addr, portnum));
  } else if(proto == "raw_udp") {
    pimpl = shared_ptr<detail::raw_udp_client_impl>(new detail::raw_udp_client_impl(ip4addr, portnum));
  };
  readNames();
  currentColumn = 0;
  reading_thread = shared_ptr<ReaKaux::thread>(new ReaKaux::thread(extract_process(*this)));
};



};


};



