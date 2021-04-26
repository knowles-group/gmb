#ifndef INIT_H
#define INIT_H

#include "get_integral.h"

#include "hamiltonian.h"
#include "supercontainer.h"

namespace bbo {
  
  void init(std::string filename, std::string method, hamiltonian<> &hamiltonian);

} // namespace bbo

#endif // INIT_H
