
#include <ReaK/core/rpc/rpc_function.hpp>

#include <iostream>

using namespace ReaK;

rpc::function< bool(int, double) > foo;


bool some_free_func(int, double) { 
  std::cout << "Hello! Free!" << std::endl;
  return true;
};

struct some_functor {
  bool operator()(int, double) const {
    std::cout << "Hello! Functor!" << std::endl;
    return false;
  };
};


int main() {
  
  // Set the 'foo' function to point to a call made to IP address 127.0.0.1
  foo.from_remote("foo", "127.0.0.1");
  
  // Call the foo function like any other function (like std::function):
//   bool b = foo(42, double());
  
  // Set the 'foo' function to point to a local function, and advertize it on a globally-set port (e.g. 17070):
  foo.publish("foo", some_free_func);
  
  // Call the foo function like any other function (like std::function), even if it is being published:
  bool b = foo(42, double());
  
  // publish should work with anything that can be converted to std::function, e.g.,:
  foo.publish("foo", 
  [](int, double) -> bool { 
    std::cout << "Hello! Lambda!" << std::endl;
    return false;
  });
  
  b = foo(42, double());
  
  foo.publish("foo", some_functor());
  
  b = foo(42, double());
  
  return 0;
};


