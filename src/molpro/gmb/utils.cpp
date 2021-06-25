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

template void get_polval(std::stringstream &ss, unsigned int &value); 
template void get_polval(std::stringstream &ss, double &value); 

} // namespace gmb

