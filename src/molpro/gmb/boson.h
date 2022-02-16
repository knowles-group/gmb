#ifndef BOSON_H
#define BOSON_H

#include <string>
#include <cmath>

using sym_t=unsigned int;
using syms_t=std::vector<sym_t>; ///< position in orbital space for each symmetry


struct polariton {
  bool self_energy{true};
  bool coupling{true};
  sym_t nmax{0};
  double gamma{0.0};
  double omega{0.0};
  double lambd{0.0};
  std::string fname_dm;
  std::string fname_sm;

  polariton() = default;
  polariton(sym_t nmax_, double gamma_, double omega_, bool self_energy_, bool coupling_) 
  : nmax{nmax_}, gamma{gamma_}, omega{omega_}, self_energy{self_energy_}, coupling{coupling_} 
  {
    lambd = gamma*(sqrt(2*omega));
  }

};

struct vibration {
  bool coupling{true};
  sym_t nmax{0};
  double omega{0.0};
  std::vector<std::string> integral_files;

  vibration() = default;
  vibration(sym_t nmax_, double omega_) 
  : nmax{nmax_}, omega{omega_}
  {
    integral_files.resize(3);
  }

};


struct bosons {
  std::vector<std::unique_ptr<polariton>> polaritons;
  std::vector<std::unique_ptr<vibration>> vibrations;
};



#endif // BOSON_H
