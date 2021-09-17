
#include "utils.h"
#include <fstream>

namespace gmb {

  /**
   * @brief Check if a file exists.
   * 
   * If file cannot be found, terminates program.
   * 
   * @param file 
   */
  void check_file(const std::string &file, const std::string &str) {
    std::ifstream infile;
    infile.open(file);
    //send error if file not found
    if (!infile) {
      if (file.size() == 0) 
        std::cerr << "No " << str << " file was provided.\n";
      else 
        std::cerr << "Unable to open " << str << " file: \"" << file << "\"\n";
      exit(1); // terminate with error
    }
    infile.close();
  }

  size_t get_offset(size_t i, size_t j, size_t k, size_t l, size_t nj, size_t nk, size_t nl) {
    size_t offset = i*(nj*nk*nl) + j*(nk*nl) + k*nl + l;
    return offset;
  }

  size_t get_offset(size_t i, size_t j, size_t nj) {
    size_t offset = i*nj + j;
    return offset;
  }

  std::string tospin(size_t spin) {
    switch (spin) {
    case 0: return "a";
      break;
    case 1: return "b";
      break;
    default: return &"p"[spin]-2;
      break;
    }
  }


} // namespace gmb

