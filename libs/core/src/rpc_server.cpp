
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

#include <ReaK/core/rpc/detail/disc_headers.hpp>
#include <ReaK/core/rpc/detail/rpc_headers.hpp>
#include <ReaK/core/rpc/rpc_exceptions.hpp>
#include <ReaK/core/rpc/rpc_server.hpp>

#include <ReaK/core/base/global_rng.hpp>
#include <ReaK/core/serialization/bin_archiver.hpp>
#include <ReaK/core/serialization/protobuf_archiver.hpp>
#include <ReaK/core/serialization/xml_archiver.hpp>
#include <memory>
#include <thread>

#include <boost/asio.hpp>
#include <boost/circular_buffer.hpp>
#include <functional>
#include <set>

#include <chrono>
#include <utility>

namespace ReaK::rpc {

#if 0
std::recursive_mutex log_mtx;
#define RK_RPC_LOG(X)                                                          \
  {                                                                            \
    std::unique_lock<std::recursive_mutex> _(log_mtx);                         \
    auto t = std::chrono::high_resolution_clock::now().time_since_epoch();     \
    std::cout                                                                  \
        << "ln: " << std::setw(5) << __LINE__ << " thd: " << std::hex          \
        << std::this_thread::get_id() << std::dec << " time: " << std::setw(4) \
        << (std::chrono::duration_cast<std::chrono::seconds>(t).count() %      \
            1000)                                                              \
        << "." << std::setfill('0') << std::setw(6)                            \
        << (std::chrono::duration_cast<std::chrono::microseconds>(t).count() % \
            1000000)                                                           \
        << std::setfill(' ') << " message: " << X << std::endl;                \
  };
#else
#define RK_RPC_LOG(X)
#endif

static volatile bool server_is_alive = true;

static std::unique_ptr<serialization::iarchive> make_iarchiver(
    std::istream& in, msg_format fmt = xml_format) {
  std::unique_ptr<serialization::iarchive> p_ar;
  switch (fmt) {
    case binary_format:
      p_ar = std::make_unique<serialization::bin_iarchive>(in);
      break;
    case xml_format:
    default:
      p_ar = std::make_unique<serialization::xml_iarchive>(in);
      break;
    case protobuf_format:
      p_ar = std::make_unique<serialization::protobuf_iarchive>(in);
      break;
  }
  return p_ar;
}

static std::unique_ptr<serialization::oarchive> make_oarchiver(
    std::ostream& out, msg_format fmt = xml_format) {
  std::unique_ptr<serialization::oarchive> p_ar;
  switch (fmt) {
    case binary_format:
      p_ar = std::make_unique<serialization::bin_oarchive>(out);
      break;
    case xml_format:
    default:
      p_ar = std::make_unique<serialization::xml_oarchive>(out);
      break;
    case protobuf_format:
      p_ar = std::make_unique<serialization::protobuf_oarchive>(out);
      break;
  }
  return p_ar;
}

static std::string gen_random_name() {
  static const char* alphanum =
      "0123456789"
      "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      "abcdefghijklmnopqrstuvwxyz";

  auto& gen_rd = get_global_rng();
  std::string result(9, '0');
  for (int i = 0; i < 9; ++i) {
    result[i] = alphanum[gen_rd() % (sizeof(alphanum) - 1)];
  }
  return result;
}

namespace {

struct sort_published_functions {
  bool operator()(detail::remote_function* lhs,
                  detail::remote_function* rhs) const {
    return ((lhs->name < rhs->name) ||
            ((lhs->name == rhs->name) &&
             (lhs->get_params_hash() < rhs->get_params_hash())));
  }
};
using pub_func_set =
    std::set<detail::remote_function*, sort_published_functions>;

struct server_config_data {
  std::string name;
  boost::asio::ip::tcp::endpoint endpoint;

  explicit server_config_data(std::string aName = "")
      : name(std::move(aName)) {}

  bool operator<(const server_config_data& rhs) const {
    return this->name < rhs.name;
  }
};
using server_set = std::set<server_config_data>;

struct compare_server_endpoints {
  bool operator()(const boost::asio::ip::tcp::endpoint& lhs,
                  const boost::asio::ip::tcp::endpoint& rhs) const {
    if (lhs.address() < rhs.address()) {
      return true;
    }
    if ((lhs.address() == rhs.address()) && (lhs.port() < rhs.port())) {
      return true;
    }
    return false;
  }
};
using server_endpt_map =
    std::map<boost::asio::ip::tcp::endpoint, server_set::iterator,
             compare_server_endpoints>;

struct pending_return {
  detail::remote_function* p_func;
  std::promise<detail::call_results> out_promise;
  detail::call_preparations pre_data;
  boost::asio::ip::tcp::socket socket;
  boost::asio::ip::tcp::endpoint endpoint;

  pending_return(boost::asio::io_service& io_service,
                 detail::remote_function* aPFunc,
                 detail::call_preparations&& aPreData)
      : p_func(aPFunc), pre_data(std::move(aPreData)), socket(io_service) {}

  pending_return(pending_return&& rhs) noexcept
      : p_func(rhs.p_func),
        out_promise(std::move(rhs.out_promise)),
        pre_data(std::move(rhs.pre_data)),
        socket(std::move(rhs.socket)),
        endpoint(std::move(rhs.endpoint)) {}

  pending_return& operator=(pending_return&& rhs) noexcept {
    p_func = rhs.p_func;
    out_promise = std::move(rhs.out_promise);
    pre_data = std::move(rhs.pre_data);
    socket = std::move(rhs.socket);
    endpoint = std::move(rhs.endpoint);
    return *this;
  }

  void handle_connect(const boost::system::error_code& e) {
    RK_RPC_LOG("Handling a connection");
    if (!e) {

      // Produce the call packet.
      boost::asio::streambuf temp_out_buf;
      std::ostream temp_out(&temp_out_buf);
      auto pre_data_dump = pre_data.ss->str();
      temp_out << detail::generate_rpc_header(detail::rpc_call, pre_data.fmt,
                                              pre_data_dump.size());
      temp_out << pre_data_dump;

      RK_RPC_LOG("Sending the call packet");
      // Send the call packet.
      temp_out_buf.consume(socket.send(temp_out_buf.data()));
      RK_RPC_LOG("Successfully sent " << snt_sz << " bytes.");

      socket.close();
    } else {
      RK_RPC_LOG("An error occurred when given a connection to handle for seq: "
                 << pre_data.call_seq << " With message: " << e.message());
      socket.async_connect(endpoint,
                           [this](const boost::system::error_code& e) {
                             this->handle_connect(e);
                           });
    }
  }

  void async_connect() {
    RK_RPC_LOG("Triggering an async connection for the call to "
               << endpoint.address().to_string() << " " << endpoint.port());
    socket.async_connect(endpoint, [this](const boost::system::error_code& e) {
      this->handle_connect(e);
    });
  }
};

}  // namespace

struct server::impl {

  msg_format preferred_fmt;
  std::thread listening_thd;
  boost::asio::io_service io_service;
  boost::asio::ip::tcp::acceptor acceptor;

  std::atomic<unsigned int> cur_seq_number;

  std::mutex pub_funcs_mtx;
  pub_func_set pub_funcs;

  std::vector<std::future<void>> future_pool;

  std::mutex pending_connects_mtx;
  std::multimap<std::string, pending_return> pending_connects;

  std::mutex pending_calls_mtx;
  std::map<std::size_t, pending_return> pending_calls;

  struct call_executer {
    impl* parent;
    msg_format fmt;
    std::string payload;
    boost::asio::ip::tcp::endpoint sender_endpt;

    call_executer(impl* aParent, msg_format aFmt, std::string&& aPayload,
                  boost::asio::ip::tcp::endpoint aSenderEndPt)
        : parent(aParent),
          fmt(aFmt),
          payload(std::move(aPayload)),
          sender_endpt(std::move(aSenderEndPt)) {}

    call_executer(call_executer&& rhs) noexcept
        : parent(rhs.parent),
          fmt(rhs.fmt),
          payload(std::move(rhs.payload)),
          sender_endpt(std::move(rhs.sender_endpt)) {}

    call_executer& operator=(call_executer&& rhs) noexcept {
      parent = rhs.parent;
      fmt = rhs.fmt;
      payload = std::move(rhs.payload);
      sender_endpt = std::move(rhs.sender_endpt);
      return *this;
    }

    void operator()() {

      RK_RPC_LOG("Starting to execute call...");

      boost::asio::ip::tcp::socket socket(parent->io_service);

      std::stringstream ss(payload);
      auto p_ari = make_iarchiver(ss, fmt);
      detail::rpc_call_header r_hdr;
      (*p_ari) >> r_hdr;

      // Find the function to be called.
      detail::dummy_remote_function candidate(r_hdr.name, r_hdr.params_hash);
      auto it = parent->pub_funcs.find(&candidate);
      if (it == parent->pub_funcs.end()) {

        RK_RPC_LOG("Could not find the function among candidates.");

        // Send a "unrecognized" message back.
        std::stringstream temp_out;
        auto p_aro = make_oarchiver(temp_out, fmt);
        detail::rpc_unrecognized_header reply_hdr(r_hdr.name, r_hdr.params_hash,
                                                  r_hdr.call_seq);
        (*p_aro) << reply_hdr;

        // Produce the return packet.
        boost::asio::streambuf out_buf;
        std::ostream out(&out_buf);
        auto res_data =
            temp_out.str();  // TODO Is there a way to avoid double copy?
        out << detail::generate_rpc_header(detail::rpc_unrecognized, fmt,
                                           res_data.size());
        out << res_data;

        RK_RPC_LOG("Connecting to sender at "
                   << sender_endpt.address().to_string() << ":"
                   << r_hdr.reply_port);

        // Send the return packet.
        boost::system::error_code ec;
        sender_endpt.port(r_hdr.reply_port);
        socket.connect(sender_endpt, ec);
        if (!ec) {
          out_buf.consume(socket.send(out_buf.data()));
        }
        RK_RPC_LOG("Sent the UNRECOGNIZED packet");
        socket.close();
        return;
      }

      RK_RPC_LOG("Called function was found");

      try {
        // Make the call.
        std::stringstream temp_out;
        auto p_aro = make_oarchiver(temp_out, fmt);
        detail::rpc_return_header reply_hdr(r_hdr.name, r_hdr.params_hash,
                                            r_hdr.call_seq);
        (*p_aro) << reply_hdr;
        RK_RPC_LOG("Wrote header into the payload");
        (*it)->execute(*p_ari, *p_aro);
        RK_RPC_LOG("Executed the function");

        // Produce the return packet.
        boost::asio::streambuf out_buf;
        std::ostream out(&out_buf);
        auto res_data =
            temp_out.str();  // TODO Is there a way to avoid double copy?
        out << detail::generate_rpc_header(detail::rpc_return, fmt,
                                           res_data.size());
        out << res_data;

        RK_RPC_LOG("Connecting to sender at "
                   << sender_endpt.address().to_string() << ":"
                   << r_hdr.reply_port);

        // Send the return packet.
        boost::system::error_code ec;
        sender_endpt.port(r_hdr.reply_port);
        socket.connect(sender_endpt, ec);
        if (!ec) {
          out_buf.consume(socket.send(out_buf.data()));
        }
        RK_RPC_LOG("Sent RETURN packet");
        socket.close();
      } catch (std::exception& e) {
        RK_RPC_LOG("Exception occurred during call");
        // Send a "exception" message back.
        std::stringstream temp_out;
        auto p_aro = make_oarchiver(temp_out, fmt);
        detail::rpc_exception_header reply_hdr(
            r_hdr.name, r_hdr.params_hash, r_hdr.call_seq,
            "Execution error on remote host", e.what());
        (*p_aro) << reply_hdr;

        // Produce the return packet.
        boost::asio::streambuf out_buf;
        std::ostream out(&out_buf);
        auto res_data =
            temp_out.str();  // TODO Is there a way to avoid double copy?
        out << detail::generate_rpc_header(detail::rpc_exception, fmt,
                                           res_data.size());
        out << res_data;

        RK_RPC_LOG("Connecting to sender at "
                   << sender_endpt.address().to_string() << ":"
                   << r_hdr.reply_port);

        // Send the return packet.
        boost::system::error_code ec;
        sender_endpt.port(r_hdr.reply_port);
        socket.connect(sender_endpt, ec);
        if (!ec) {
          out_buf.consume(socket.send(out_buf.data()));
        }
        RK_RPC_LOG("Sent the EXCEPTION packet");
        socket.close();
      }
    }
  };

  struct connection {
    impl* parent;

    boost::asio::ip::tcp::iostream socket_stream;
    boost::asio::ip::tcp::endpoint endpoint;

    explicit connection(impl* aParent) : parent(aParent){};

    bool check_call_seq_match(
        std::map<std::size_t, pending_return>::iterator it,
        const detail::rpc_basic_header& r_hdr) const {
      if ((it->second.p_func->name != r_hdr.name) ||
          (it->second.p_func->get_params_hash() != r_hdr.params_hash)) {
        std::stringstream ss;
        ss << "The remote call " << r_hdr.call_seq
           << " returned a mismatch of function! Expected '"
           << it->second.p_func->name << ":"
           << it->second.p_func->get_params_hash() << "' but got '"
           << r_hdr.name << ":" << r_hdr.params_hash << "'!" << std::flush;
        it->second.out_promise.set_exception(
            std::make_exception_ptr(communication_error(ss.str())));
        parent->pending_calls.erase(it);
        return false;
      };
      return true;
    };

    void operator()(std::shared_ptr<connection> /* to keep self alive */,
                    const boost::system::error_code& e) {

      parent->start_accept();

      {
        if (e.value() != 0) {
          RK_RPC_LOG("Error on tcp socket: '"
                     << e.message()
                     << "' with error-code value = " << e.value());
        } else {
          RK_RPC_LOG("Connection on tcp socket was empty!");
        }
      }

      detail::rpc_header_type hdr_type = {};
      msg_format fmt = {};
      std::size_t msg_size = 0;

      while (socket_stream) {
        try {

          std::tie(hdr_type, fmt, msg_size) =
              detail::parse_rpc_header(socket_stream);
          RK_RPC_LOG("Received RPC packet: "
                     << detail::rpc_header_type_to_str[hdr_type] << " of size "
                     << msg_size);

          // Check for either request for a procedure execution
          //  or reply from a remotely executed call.
          switch (hdr_type) {
            case detail::rpc_call: {

              auto it = parent->future_pool.begin();
              for (; it != parent->future_pool.end(); ++it) {
                if (it->wait_for(std::chrono::nanoseconds(0)) ==
                    std::future_status::ready) {
                  break;
                }
              }

              // Defer rpc call to a separate thread:
              std::string payload(msg_size, '\0');
              socket_stream.read(payload.data(), msg_size);
              auto f = std::async(
                  std::launch::async,
                  call_executer(parent, fmt, std::move(payload), endpoint));
              if (it == parent->future_pool.end()) {
                parent->future_pool.push_back(std::move(f));
              } else {
                *it = std::move(f);
              }

              break;
            }
            case detail::rpc_return: {
              // Must copy data from socket_stream before it dies.
              detail::call_results res;
              std::string payload(
                  msg_size, '\0');  // TODO Is there a way to avoid double copy?
              socket_stream.read(payload.data(), msg_size);
              res.ss = std::make_unique<std::stringstream>(payload);
              res.p_ari = make_iarchiver(*res.ss, fmt);
              detail::rpc_return_header r_hdr;
              (*res.p_ari) >> r_hdr;

              // Find the function-call promise in question.
              std::unique_lock<std::mutex> _(parent->pending_calls_mtx);
              auto it = parent->pending_calls.find(r_hdr.call_seq);
              if (!check_call_seq_match(it, r_hdr)) {
                break;
              }

              // Fulfill the function-call's promise by delivering the input-archive.
              it->second.out_promise.set_value(std::move(res));
              parent->pending_calls.erase(it);
              break;
            }
            case detail::rpc_exception: {
              auto p_ari = make_iarchiver(socket_stream, fmt);
              detail::rpc_exception_header r_hdr;
              (*p_ari) >> r_hdr;

              // Find the function-call promise in question.
              std::unique_lock<std::mutex> _(parent->pending_calls_mtx);
              auto it = parent->pending_calls.find(r_hdr.call_seq);
              if (!check_call_seq_match(it, r_hdr)) {
                break;
              }

              // Translate the rpc-exception to regular exception.
              auto e_p = std::make_exception_ptr(std::runtime_error(
                  "RPC Server reported an exception of type '" +
                  r_hdr.except_type + "' with message: '" + r_hdr.except_msg +
                  "'."));
              /* TODO Maybe translate the exception better */

              // Issue the exception on the function-call's promise.
              it->second.out_promise.set_exception(e_p);
              parent->pending_calls.erase(it);
              break;
            }
            case detail::rpc_unrecognized: {
              auto p_ari = make_iarchiver(socket_stream, fmt);
              detail::rpc_unrecognized_header r_hdr;
              (*p_ari) >> r_hdr;

              // Find the function-call promise in question.
              std::unique_lock<std::mutex> _(parent->pending_calls_mtx);
              auto it = parent->pending_calls.find(r_hdr.call_seq);
              if (!check_call_seq_match(it, r_hdr)) {
                break;
              }

              // Issue the unrecognized exception on the function-call's promise.
              it->second.out_promise.set_exception(std::make_exception_ptr(
                  unrecognized_function(it->second.p_func)));
              parent->pending_calls.erase(it);
              break;
            }
            case detail::rpc_startportserver: {
              parent->discovery_hdl.start_port_server();
              break;
            }
            default:
              break;
          }

        } catch (...) {
          // nothing to do if this fails (bad packet, bad whatever.. just can't handle the packet, so, leave it)
          /* TODO Log the error (maybe) */
        }
      }
    }
  };

  struct discovery_handler {
   private:
    impl* parent;

    boost::asio::ip::udp::socket udp_socket;
    boost::asio::ip::udp::endpoint multicast_endpt;
    boost::asio::ip::udp::endpoint sender_endpt;
    boost::asio::streambuf buf;

    std::mutex self_config_mtx;
    std::atomic<bool> self_is_port_server;
    server_config_data self_config;
    server_set all_servers;
    server_endpt_map assigned_endpts;
    std::atomic<bool> self_config_confirmed;

    std::promise<bool> given_addr_port;
    std::future<void> startup_executer;

   public:
    explicit discovery_handler(impl* aParent)
        : parent(aParent),
          udp_socket(parent->io_service),
          multicast_endpt(boost::asio::ip::address::from_string("239.255.17.1"),
                          17001),
          self_is_port_server(false),
          self_config(gen_random_name()),
          self_config_confirmed(false) {

      // Temporarily set local endpoint to whatever:
      self_config.endpoint =
          boost::asio::ip::tcp::endpoint(boost::asio::ip::tcp::v4(), 17000);

      // Create a socket for receiving udp multicasting messages from discovery protocol.
      boost::asio::ip::udp::endpoint listen_endpoint(boost::asio::ip::udp::v4(),
                                                     17001);
      udp_socket.open(listen_endpoint.protocol());
      udp_socket.set_option(boost::asio::ip::udp::socket::reuse_address(true));
      udp_socket.set_option(boost::asio::ip::multicast::enable_loopback(true));
      udp_socket.bind(listen_endpoint);

      // Join the multicast group.
      udp_socket.set_option(
          boost::asio::ip::multicast::join_group(multicast_endpt.address()));

      udp_socket.async_receive_from(
          buf.prepare(1024), sender_endpt,
          [this](const boost::system::error_code& e, std::size_t bytes_recvd) {
            handle_receive_from(e, bytes_recvd);
          });

      startup_executer = std::async(std::launch::async, [this]() {
        startup_procedure(self_config.name,
                          self_config.endpoint.address().to_string());
      });
    }

    void finalize() {
      // NOTE finalize is called on the io_service thread (after stop)
      //      this means that only synchronous operations can be done here
      //      but it also means that locks could be elided where appropriate.

      // first, I need to make sure self is active (has a valid port assigned to it):
      if (self_config.endpoint.port() > 17000) {
        std::unique_lock<std::mutex> _(self_config_mtx);
        // IF (self is PSv) THEN
        if (self_is_port_server) {
          boost::asio::ip::tcp::socket socket(parent->io_service);
          socket.set_option(boost::asio::ip::tcp::socket::reuse_address(true));
          boost::system::error_code asio_err =
              boost::asio::error::host_not_found;
          // mark self as NOT PSv
          self_is_port_server = false;
          // select any active server
          for (const auto& sv : all_servers) {  // DO
            if (sv.name == self_config.name) {
              continue;
            }

            // send STARTPORTSERVER to selected server
            try {
              // Produce the startportserver packet.
              boost::asio::streambuf temp_out_buf;
              std::ostream temp_out(&temp_out_buf);
              temp_out << detail::generate_rpc_header(
                  detail::rpc_startportserver, parent->preferred_fmt, 0);

              // Prepare the socket.
              socket.close();
              socket.connect(sv.endpoint, asio_err);
              if (asio_err) {
                continue;
              }

              // Send the call packet.
              temp_out_buf.consume(socket.send(temp_out_buf.data()));
              break;
            } catch (std::exception& e) {
              RK_UNUSED(e);
              continue;
            }
          }  // UNTIL (TCP packet was sent)
        }    // END

        // broadcast UNADVERTISE
        unadvertise_bad_port();

        // remove self addr-port
        self_config.endpoint.port(17000);

        // stop TCP listening on addr+port
        rebind_acceptor();
      }

      udp_socket.cancel();

      try {
        given_addr_port.set_value(true);
      } catch (...) {}

      startup_executer.get();

      udp_socket.close();
    }

    void start_port_server() { self_is_port_server.store(true); };

    void set_self_name(const std::string& aName) {
      std::unique_lock<std::mutex> _(self_config_mtx);

      self_config.name = aName;

      if (self_config.endpoint.port() > 17001) {
        advertise_new_port();
      }
    }

    void rebind_acceptor() {
      if (parent->acceptor.is_open()) {
        parent->acceptor.cancel();
        parent->acceptor.close();
      }
      if (self_config.endpoint.port() > 17000) {
        parent->acceptor = boost::asio::ip::tcp::acceptor(
            parent->io_service,
            boost::asio::ip::tcp::endpoint(boost::asio::ip::tcp::v4(),
                                           self_config.endpoint.port()));
      }
    }

    void wait_for_discovery() const {
      // Spin-lock on the validity of the port:
      while (!self_config_confirmed.load()) {
        std::this_thread::yield();
      }
    }

    unsigned short get_self_port() {
      std::unique_lock<std::mutex> _(self_config_mtx);
      return self_config.endpoint.port();
    }

    std::future<detail::call_results> make_remote_call(
        pending_return&& pre_connection) {
      bool got_ep = false;
      {
        server_config_data sv_data(pre_connection.p_func->get_host());
        std::unique_lock<std::mutex> _(self_config_mtx);
        RK_RPC_LOG("Took the lock, now looking for requested host");
        auto it = all_servers.find(sv_data);
        if (it != all_servers.end()) {
          RK_RPC_LOG("Found the host! at address: "
                     << it->endpoint.address().to_string() << " "
                     << it->endpoint.port());
          pre_connection.endpoint = it->endpoint;
          if (pre_connection.endpoint.address() ==
              self_config.endpoint.address()) {
            pre_connection.endpoint.address(
                boost::asio::ip::address::from_string("127.0.0.1"));
          }
          got_ep = true;
        }
      }
      if (got_ep) {
        // Add to the pending calls
        std::unique_lock<std::mutex> _(parent->pending_calls_mtx);
        RK_RPC_LOG("Took lock on pending calls queue");
        auto sq = pre_connection.pre_data.call_seq;
        auto it = parent->pending_calls
                      .insert(std::make_pair(sq, std::move(pre_connection)))
                      .first;
        it->second.async_connect();
        RK_RPC_LOG(
            "Triggered the async connection for the pending call, returning "
            "future..");
        return it->second.out_promise.get_future();
      }
      std::unique_lock<std::mutex> _(parent->pending_connects_mtx);
      RK_RPC_LOG("Took a lock on connection queue");
      auto* pf = pre_connection.p_func;
      auto it = parent->pending_connects.insert(
          std::make_pair(pf->get_host(), std::move(pre_connection)));
      RK_RPC_LOG(
          "Added pre-connection object to the connection queue, returning "
          "future..");
      return it->second.out_promise.get_future();
    }

   private:
    void startup_procedure(const std::string& aSelfName,
                           const std::string& aSelfAddr) {
      try {

        // Make the NEWPORT broadcast call.
        std::stringstream temp_out;
        auto p_aro = make_oarchiver(temp_out, parent->preferred_fmt);
        detail::disc_basic_header np_hdr(aSelfName, aSelfAddr);
        (*p_aro) << np_hdr;

        // Produce the newport request packet.
        boost::asio::streambuf out_buf;
        std::ostream out(&out_buf);
        std::string res_data = temp_out.str();
        out << detail::generate_disc_header(
            detail::disc_newport, parent->preferred_fmt, res_data.size());
        out << res_data;

        auto obtained_port = given_addr_port.get_future();

        // Send the newport request packet.
        out_buf.consume(udp_socket.send_to(out_buf.data(), multicast_endpt));

        if (obtained_port.wait_for(std::chrono::seconds(1)) !=
            std::future_status::ready) {
          // claim the role of the port server:
          self_is_port_server.store(true);
        }

      } catch (std::exception& e) {
        RK_UNUSED(e);
        // TODO Do something, what? FIXME
        RK_RPC_LOG(
            "Failed the send the NEWPORT request... exception: " << e.what());
      }
    }

    void delayed_startup_procedure(const std::string& aSelfName,
                                   const std::string& aSelfAddr) {
      std::this_thread::sleep_for(std::chrono::seconds(1));
      startup_procedure(aSelfName, aSelfAddr);
    }

    void advertise_new_port() {
      // Make the ADVERTISE broadcast call.
      std::stringstream temp_out;
      auto p_aro = make_oarchiver(temp_out, parent->preferred_fmt);
      detail::disc_with_port_header wp_hdr(
          self_config.name, self_config.endpoint.address().to_string(),
          self_config.endpoint.port());
      (*p_aro) << wp_hdr;

      // Produce the advertise request packet.
      boost::asio::streambuf out_buf;
      std::ostream out(&out_buf);
      std::string res_data = temp_out.str();
      out << detail::generate_disc_header(
          detail::disc_advertise, parent->preferred_fmt, res_data.size());
      out << res_data;

      // Send the advertise request packet.
      out_buf.consume(udp_socket.send_to(out_buf.data(), multicast_endpt));
    }

    void echo_self_port() {
      // Make the ECHO broadcast call.
      std::stringstream temp_out;
      auto p_aro = make_oarchiver(temp_out, parent->preferred_fmt);
      detail::disc_with_port_header wp_hdr(
          self_config.name, self_config.endpoint.address().to_string(),
          self_config.endpoint.port());
      (*p_aro) << wp_hdr;

      // Produce the echo request packet.
      boost::asio::streambuf out_buf;
      std::ostream out(&out_buf);
      std::string res_data = temp_out.str();
      out << detail::generate_disc_header(
          detail::disc_echo, parent->preferred_fmt, res_data.size());
      out << res_data;

      // Send the echo request packet.
      out_buf.consume(udp_socket.send_to(out_buf.data(), multicast_endpt));
    }

    void unadvertise_bad_port() {
      if (self_is_port_server) {
        self_is_port_server = false;
      }
      // Make the UNADVERTISE broadcast.
      std::stringstream temp_out;
      auto p_aro = make_oarchiver(temp_out, parent->preferred_fmt);
      detail::disc_basic_header np_hdr(
          self_config.name, self_config.endpoint.address().to_string());
      (*p_aro) << np_hdr;

      // Produce the unadvertise request packet.
      boost::asio::streambuf out_buf;
      std::ostream out(&out_buf);
      std::string res_data = temp_out.str();
      out << detail::generate_disc_header(
          detail::disc_unadvertise, parent->preferred_fmt, res_data.size());
      out << res_data;

      // Send the unadvertise request packet.
      out_buf.consume(udp_socket.send_to(out_buf.data(), multicast_endpt));
    }

    void assign_addr_port_advertised(
        const detail::disc_with_port_header& r_hdr) {
      server_config_data tmp_sv_data(r_hdr.name);
      tmp_sv_data.endpoint.address(
          boost::asio::ip::address::from_string(r_hdr.addr));
      tmp_sv_data.endpoint.port(static_cast<unsigned short>(r_hdr.port));

      {
        std::unique_lock<std::mutex> _(self_config_mtx);
        auto old_ep_it = assigned_endpts.find(tmp_sv_data.endpoint);
        if (old_ep_it != assigned_endpts.end()) {
          all_servers.erase(old_ep_it->second);
          assigned_endpts.erase(old_ep_it);
        }
        auto it = all_servers.find(tmp_sv_data);
        if (it != all_servers.end()) {
          auto ep_it = assigned_endpts.find(it->endpoint);
          if (ep_it != assigned_endpts.end()) {
            assigned_endpts.erase(ep_it);
          }
          all_servers.erase(it);
        }
        std::tie(it, std::ignore) = all_servers.insert(tmp_sv_data);
        assigned_endpts[tmp_sv_data.endpoint] = it;
      }

      if (tmp_sv_data.endpoint.address() == self_config.endpoint.address()) {
        tmp_sv_data.endpoint.address(
            boost::asio::ip::address::from_string("127.0.0.1"));
      }

      // See if the newly advertised server is something that was pending a call:
      std::unique_lock<std::mutex> _(parent->pending_connects_mtx);
      auto rg = parent->pending_connects.equal_range(tmp_sv_data.name);
      if (rg.first != rg.second) {
        // NOTE I feel a little uneasy about holding 2 locks simultaneously.
        std::unique_lock<std::mutex> _l1(parent->pending_calls_mtx);
        while (rg.first != rg.second) {
          rg.first->second.endpoint = tmp_sv_data.endpoint;
          // Add to the pending calls
          auto sq = rg.first->second.pre_data.call_seq;
          auto it2 =
              parent->pending_calls
                  .insert(std::make_pair(sq, std::move(rg.first->second)))
                  .first;
          it2->second.async_connect();
          rg.first = parent->pending_connects.erase(rg.first);
        }
      }
    }

    void remove_name_unadvertised(const detail::disc_basic_header& r_hdr) {
      server_config_data tmp_sv_data(r_hdr.name);
      auto it = all_servers.find(tmp_sv_data);
      if (it != all_servers.end()) {
        auto ep_it = assigned_endpts.find(it->endpoint);
        if (ep_it != assigned_endpts.end()) {
          assigned_endpts.erase(ep_it);
        }
        all_servers.erase(it);
      }
    }

    void handle_receive_from(const boost::system::error_code& e,
                             std::size_t bytes_recvd) {
      if (!e && (bytes_recvd != 0)) {
        try {

          buf.commit(bytes_recvd);

          detail::disc_header_type hdr_type = {};
          msg_format fmt = {};
          std::size_t msg_size = 0;
          std::istream in(&buf);

          std::tie(hdr_type, fmt, msg_size) = detail::parse_disc_header(in);

          auto p_ari = make_iarchiver(in, fmt);

          switch (hdr_type) {
            case detail::disc_newport: {
              detail::disc_basic_header r_hdr;
              (*p_ari) >> r_hdr;
              RK_RPC_LOG("Got NEWPORT from " << r_hdr.name << " with addr "
                                             << r_hdr.addr);

              if (r_hdr.name == self_config.name) {
                if (sender_endpt.address() != self_config.endpoint.address()) {
                  // detected a loopback of one's own newport request, use that to get local IP address:
                  RK_RPC_LOG("Self-recognized address as "
                             << sender_endpt.address());
                  std::unique_lock<std::mutex> _(self_config_mtx);
                  self_config.endpoint.address(sender_endpt.address());
                }
              } else if (self_is_port_server) {
                std::unique_lock<std::mutex> _(self_config_mtx);

                // find free port -> N
                unsigned short new_port = 17002;
                bool must_advertise_self = false;
                if (self_config.endpoint.port() == 17000) {
                  RK_RPC_LOG("This is the first server to join network");
                  self_config.endpoint.port(17002);
                  new_port = 17003;
                  must_advertise_self = true;
                } else {
                  RK_RPC_LOG("Assigning port to the server");

                  server_config_data tmp_sv_data(r_hdr.name);
                  auto it = all_servers.find(tmp_sv_data);
                  if (it != all_servers.end()) {
                    break;  // already has a port.
                  }

                  RK_RPC_LOG("Looking for free port for the server");

                  auto ep_lo = assigned_endpts.lower_bound(
                      boost::asio::ip::tcp::endpoint(sender_endpt.address(),
                                                     17002));
                  auto ep_up = assigned_endpts.upper_bound(
                      boost::asio::ip::tcp::endpoint(sender_endpt.address(),
                                                     17999));
                  while (ep_lo != ep_up) {
                    if (ep_lo->first.port() > new_port) {
                      break;
                    }
                    ++ep_lo;
                    ++new_port;
                  }
                }

                RK_RPC_LOG("Assigned the port " << new_port
                                                << " for the server");

                // send GIVEPORT "N" to source of request
                std::stringstream temp_out;
                auto p_aro = make_oarchiver(temp_out, parent->preferred_fmt);
                detail::disc_with_port_header wp_hdr(
                    r_hdr.name, sender_endpt.address().to_string(), new_port);
                (*p_aro) << wp_hdr;

                // Produce the giveport packet.
                boost::asio::streambuf out_buf;
                std::ostream out(&out_buf);
                std::string res_data = temp_out.str();
                out << detail::generate_disc_header(detail::disc_giveport,
                                                    parent->preferred_fmt,
                                                    res_data.size());
                out << res_data;

                // Send the giveport packet.
                out_buf.consume(
                    udp_socket.send_to(out_buf.data(), multicast_endpt));

                RK_RPC_LOG("Sent the GIVEPORT packet to the server");

                // broadcast ADVERTISE
                if (must_advertise_self) {

                  // start TCP listening on addr+port
                  rebind_acceptor();
                  parent->start_accept();

                  RK_RPC_LOG("Advertising the self-assigned new port "
                             << self_config.endpoint.port());
                  // broadcast ADVERTISE
                  advertise_new_port();
                  self_config_confirmed.store(true);
                }

              } else if (self_config.endpoint.port() == 17000) {
                try {
                  throw communication_error("Concurrent NEWPORT request!");
                } catch (...) {
                  given_addr_port.set_exception(std::current_exception());
                  given_addr_port = std::promise<bool>();
                  startup_executer = std::async(std::launch::async, [this]() {
                    delayed_startup_procedure(
                        self_config.name,
                        self_config.endpoint.address().to_string());
                  });
                }
              }

              break;
            }
            case detail::disc_giveport: {
              detail::disc_with_port_header r_hdr;
              (*p_ari) >> r_hdr;
              RK_RPC_LOG("Got GIVEPORT from " << r_hdr.name << " with addr "
                                              << r_hdr.addr << ":"
                                              << r_hdr.port);

              if (r_hdr.name != self_config.name) {
                break;
              }

              std::unique_lock<std::mutex> _(self_config_mtx);

              // IF (self already had port) THEN
              if (self_config.endpoint.port() > 17001) {
                RK_RPC_LOG("Given port overrides previously given port");

                // broadcast UNADVERTISE
                unadvertise_bad_port();

                // remove self addr+port
                self_config.endpoint.port(17000);

                // stop TCP listening on addr+port
                rebind_acceptor();
              };  // END

              RK_RPC_LOG("Adopting the given port");

              // assign self addr+port
              self_config.endpoint.port(r_hdr.port);
              self_config.endpoint.address(
                  boost::asio::ip::address::from_string(r_hdr.addr));
              given_addr_port.set_value(true);

              // start TCP listening on addr+port
              rebind_acceptor();
              parent->start_accept();

              RK_RPC_LOG("Advertising the new port");
              // broadcast ADVERTISE
              advertise_new_port();
              self_config_confirmed.store(true);

              break;
            }
            case detail::disc_advertise: {
              detail::disc_with_port_header r_hdr;
              (*p_ari) >> r_hdr;
              RK_RPC_LOG("Got ADVERTISE from " << r_hdr.name << " with addr "
                                               << r_hdr.addr << ":"
                                               << r_hdr.port);

              if ((r_hdr.name == self_config.name) &&
                  (sender_endpt.address() != self_config.endpoint.address())) {
                RK_RPC_LOG("Self-recognized address as "
                           << sender_endpt.address());
                // detected a loopback of one's own advertise request, use that to get local IP address:
                std::unique_lock<std::mutex> _(self_config_mtx);
                self_config.endpoint.address(sender_endpt.address());
              }

              // assign addr+port to advertised name
              assign_addr_port_advertised(r_hdr);

              // broadcast ECHO for self addr+port
              if ((r_hdr.name != self_config.name) &&
                  (self_config.endpoint.port() > 17001)) {
                RK_RPC_LOG("Echoing own addr-port");
                std::unique_lock<std::mutex> _(self_config_mtx);
                echo_self_port();
              }

              break;
            }
            case detail::disc_unadvertise: {
              detail::disc_basic_header r_hdr;
              (*p_ari) >> r_hdr;
              RK_RPC_LOG("Got UNADVERTISE from " << r_hdr.name << " with addr "
                                                 << r_hdr.addr);

              // remove unadvertised name + addr+port
              std::unique_lock<std::mutex> _(self_config_mtx);
              remove_name_unadvertised(r_hdr);

              break;
            }
            case detail::disc_echo: {
              detail::disc_with_port_header r_hdr;
              (*p_ari) >> r_hdr;
              RK_RPC_LOG("Got ECHO from " << r_hdr.name << " with addr "
                                          << r_hdr.addr << ":" << r_hdr.port);

              if (r_hdr.name == self_config.name) {
                break;
              }

              // assign addr+port to advertised name
              assign_addr_port_advertised(r_hdr);

              // IF (addr+port conflicts with self) THEN
              if (self_config.endpoint ==
                  boost::asio::ip::tcp::endpoint(
                      boost::asio::ip::address::from_string(r_hdr.addr),
                      r_hdr.port)) {
                RK_RPC_LOG(
                    "Found conflict with self addr-port... trying to take a "
                    "lock on self-config");
                std::unique_lock<std::mutex> _(self_config_mtx);
                // broadcast UNADVERTISE
                unadvertise_bad_port();

                // remove self addr-port
                self_config.endpoint.port(17000);

                // stop TCP listening on addr+port
                rebind_acceptor();

                // goto startup procedure
                try {
                  throw communication_error("Conflicting ECHO broadcast!");
                } catch (...) {
                  given_addr_port.set_exception(std::current_exception());
                  given_addr_port = std::promise<bool>();
                  startup_executer = std::async(std::launch::async, [this]() {
                    startup_procedure(
                        self_config.name,
                        self_config.endpoint.address().to_string());
                  });
                }
              }  // END
              break;
            }
            default:
              break;
          }

        } catch (...) {
          // nothing to do if this fails (bad packet, bad whatever.. just can't handle the packet, so, leave it)
          /* TODO Log the error (maybe) */
        }
      }

      udp_socket.async_receive_from(
          buf.prepare(1024), sender_endpt,
          [this](const boost::system::error_code& e, std::size_t bytes_recvd) {
            handle_receive_from(e, bytes_recvd);
          });
      RK_RPC_LOG("Async reception of udp discovery packets restarted");
    }
  };

  discovery_handler discovery_hdl;

  detail::call_preparations prepare_call(
      const detail::remote_function* aPFunc) {
    // Spin-lock on discovery:
    discovery_hdl.wait_for_discovery();
    // Prepare archive of parameters
    detail::call_preparations result;
    result.fmt = preferred_fmt;
    result.ss = std::make_unique<std::stringstream>();
    result.p_aro = make_oarchiver(*result.ss, result.fmt);
    result.call_seq = cur_seq_number++;
    detail::rpc_call_header call_hdr(aPFunc->name, aPFunc->get_params_hash(),
                                     result.call_seq,
                                     discovery_hdl.get_self_port());
    (*result.p_aro) << call_hdr;
    return result;
  }

  std::future<detail::call_results> make_remote_call(
      detail::remote_function* aPFunc, detail::call_preparations&& pre_data) {
    RK_RPC_LOG("Initiating remote call with seq number: " << pre_data.call_seq);
    return discovery_hdl.make_remote_call(
        pending_return(io_service, aPFunc, std::move(pre_data)));
  }

  void start_accept() {
    // NOTE I wish I could use unique_ptr instead, but no C++14 lambdas,
    //      and Boost.Asio is a retarded antiquated fucking bullshit library!!
    auto new_conn = std::make_shared<connection>(this);
    acceptor.async_accept(*(new_conn->socket_stream.rdbuf()),
                          new_conn->endpoint,
                          [new_conn](const boost::system::error_code& e) {
                            (*new_conn)(new_conn, e);
                          });
  }

  void set_self_name(const std::string& aName) {
    discovery_hdl.set_self_name(aName);
  }

  impl() : acceptor(io_service), cur_seq_number(1), discovery_hdl(this) {

    listening_thd = std::thread([this]() {
      this->io_service.run();
      this->discovery_hdl.finalize();
    });
  }

  ~impl() {
    io_service.stop();
    listening_thd.join();
  }
};

static bool are_same_remotes(detail::remote_function* lhs,
                             detail::remote_function* rhs) {
  return ((lhs->name == rhs->name) &&
          (lhs->get_params_hash() == rhs->get_params_hash()));
}

server::server() : pimpl(new server::impl()) {
  server_is_alive = true;
  // NOTE: technically, this is undefined behavior, since the 'server_is_alive' might not
  //       have been initialized yet, but since it is static (storage) and trivial type
  //       this is safe.
}

server::~server() {
  server_is_alive = false;
  delete pimpl;
}

server& server::instance() {
  static server inst;
  return inst;
}

void server::set_name(const std::string& aName) {
  instance().pimpl->set_self_name(aName);
}

void server::publish_function(detail::remote_function* aFunc) {
  if (!server_is_alive) {
    return;
  }
  std::unique_lock<std::mutex> _(pimpl->pub_funcs_mtx);
  auto it = pimpl->pub_funcs.lower_bound(aFunc);
  if ((it != pimpl->pub_funcs.end()) && are_same_remotes(*it, aFunc)) {
    if (aFunc != *it) {
      throw publishing_mismatch("publish", aFunc, *it);
    }
    return;
  }
  pimpl->pub_funcs.insert(it, aFunc);
}

void server::unpublish_function(detail::remote_function* aFunc) {
  if (!server_is_alive) {
    return;
  }
  std::unique_lock<std::mutex> _(pimpl->pub_funcs_mtx);
  auto it = pimpl->pub_funcs.lower_bound(aFunc);
  if ((it != pimpl->pub_funcs.end()) && are_same_remotes(*it, aFunc)) {
    if (aFunc != *it) {
      throw publishing_mismatch("unpublish", aFunc, *it);
    }
    pimpl->pub_funcs.erase(it);
  }
}

detail::call_preparations server::prepare_call(
    detail::remote_function* aPFunc) const {
  assert(
      server_is_alive);  // NOTE: this will only occur during static destruction, can't throw here.
  return pimpl->prepare_call(aPFunc);
}

std::future<detail::call_results> server::make_remote_call(
    detail::remote_function* aPFunc,
    detail::call_preparations&& pre_data) const {
  assert(
      server_is_alive);  // NOTE: this will only occur during static destruction, can't throw here.
  return pimpl->make_remote_call(aPFunc, std::move(pre_data));
}

}  // namespace ReaK::rpc
