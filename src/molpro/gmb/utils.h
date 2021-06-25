#ifndef GMB_UTILS_H
#define GMB_UTILS_H

#include "libtensor/libtensor.h"

// support libtensor functions 
namespace gmb {

  template <typename T>
  void get_polval(std::stringstream &ss, T &value);

  void check_file(std::string file);
    
} // namespace gmb
#endif //GMB_UTILS_H
