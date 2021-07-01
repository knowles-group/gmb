#ifndef GMB_UTILS_H
#define GMB_UTILS_H

#include "libtensor/libtensor.h"

// support libtensor functions 
namespace gmb {

  template <typename T>
  void get_polval(std::stringstream &ss, T &value);

  void check_file(std::string file);
  
  size_t get_offset(size_t i, size_t j, size_t k, size_t l, size_t nj, size_t nk, size_t nl);
  size_t get_offset(size_t i, size_t j, size_t nj);
    
} // namespace gmb
#endif //GMB_UTILS_H
