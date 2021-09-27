#ifndef GMB_SRC_MOLPRO_GMB_UTILS_H
#define GMB_SRC_MOLPRO_GMB_UTILS_H

#include <libtensor/libtensor.h>

// support libtensor functions 
namespace gmb {


  void check_file(const std::string &file, const std::string &str);
  
  size_t get_offset(size_t i, size_t j, size_t k, size_t l, size_t nj, size_t nk, size_t nl);
  size_t get_offset(size_t i, size_t j, size_t nj);

  std::string tospin(size_t spin);
    
} // namespace gmb
#endif // GMB_SRC_MOLPRO_GMB_UTILS_H
