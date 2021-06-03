#ifndef GMB_INIT_H
#define GMB_INIT_H

#include "get_integral.h"
// #include "get_integral_ph.h"
#include "get_integral_pol.h"

#include "hamiltonian.h"
#include "supercontainer.h"

namespace gmb {
  
  void init(std::string filename, std::string method, hamiltonian<> &hamiltonian);
  void init_ccpol(std::string filename, std::string method, hamiltonian<> &hamiltonian);
  // void init_ph(std::string filename, std::string method, hamiltonian<> &hamiltonian);

} // namespace gmb

#endif //GMB_INIT_H
