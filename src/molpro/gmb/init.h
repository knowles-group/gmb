#ifndef GMB_INIT_H
#define GMB_INIT_H

#include "get_integral.h"

#include "hamiltonian.h"
#include "supercontainer.h"

namespace gmb {
  
  void init(std::string filename, std::string method, hamiltonian<> &hamiltonian);

} // namespace gmb

#endif //GMB_INIT_H
