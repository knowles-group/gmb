#include "utils.h"
#include <fstream>

namespace gmb {

  template <typename T>
  void get_polval(std::stringstream &ss, T &value) {
    if (!ss.good()) {
      std::cerr << "The keyword polariton needs 3 arguments as: polariton=[nmax,gamma,omega]\n";
      exit(1);
    }
    std::string substr;
    getline(ss, substr, ',');
    std::stringstream ssval{substr};
    ssval >> value;
  }

  /**
   * @brief Check if a file exists.
   * 
   * If file cannot be found, terminates program.
   * 
   * @param file 
   */
  void check_file(std::string file) {
    std::ifstream infile;
    infile.open(file);
    //send error if file not found
    if (!infile) {
        std::cerr << "Unable to open file: " << file << "\n";
        exit(1); // terminate with error
    }
    infile.close();
  }

  size_t get_offset(size_t i, size_t j, size_t k, size_t l, size_t ni, size_t nj, size_t nk, size_t nl) {
    size_t offset = i*(nj*nk*nl) + j*(nk*nl) + k*nl + l;
    return offset;
  }

  size_t get_offset(size_t i, size_t j, size_t nj) {
    size_t offset = i*nj + j;
    return offset;
  }

template void get_polval(std::stringstream &ss, unsigned int &value); 
template void get_polval(std::stringstream &ss, double &value); 

} // namespace gmb

