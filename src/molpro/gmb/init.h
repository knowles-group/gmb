#ifndef GMB_INIT_H
#define GMB_INIT_H

#include "get_integral.h"

#include "hamiltonian.h"
#include "supercontainer.h"

namespace bbo {
  
  void init(std::string filename, std::string method, hamiltonian<> &hamiltonian);

} // namespace bbo

#endif //GMB_INIT_H
