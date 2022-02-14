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
  sym_t nmax{0};
  double omega{0.0};
  std::string fname_fock;
  // std::string fname_sm;

  vibration() = default;
  vibration(sym_t nmax_, double omega_) 
  : nmax{nmax_}, omega{omega_}
  {}

};

#endif // BOSON_H
