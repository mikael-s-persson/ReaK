
#include "ReaK/core/rpc/rpc_server.h"
#include "ReaK/core/rpc/rpc_function.h"

#include <iostream>

using namespace ReaK;

rpc::function<bool(int, double)> foo;

bool some_free_func(int /*unused*/, double /*unused*/) {
  std::cout << "Hello! Free!" << std::endl;
  return true;
};

struct some_functor {
  bool operator()(int /*unused*/, double /*unused*/) const {
    std::cout << "Hello! Functor!" << std::endl;
    return false;
  };
};

int main(int argc, char** argv) {

  if (argc < 3) {
    return 1;
  }

  std::string sv_name = argv[1];

  rpc::server::set_name(sv_name);

  rpc::function<std::string(std::string)> print_sv_name;
  rpc::function<std::string(std::string)> print_sv_name_remote;

  print_sv_name.publish("print_sv_name",
                        [&sv_name](std::string s) -> std::string {
                          std::cout << "Received call from: " << s << std::endl;
                          return sv_name;
                        });

  // Set the 'print_sv_name' function to point to a call made to host argv[2]
  print_sv_name_remote.from_remote("print_sv_name", argv[2]);

  std::cout << "Got server name from remote host as: "
            << print_sv_name_remote(sv_name) << std::endl;

  while (true) {};

  return 0;
};
