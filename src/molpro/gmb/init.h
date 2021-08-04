#ifndef GMB_INIT_H
#define GMB_INIT_H

#include "get_integral.h"
#include "utils.h"

#include "hamiltonian.h"
#include "supercontainer.h"

#include <memory>

namespace gmb {

  /**
   * @brief Initialise Hamiltonian 
   * 
   * @param filename 
   * @param method 
   * @param hamiltonian 
   */
  void init(const std::string &filename, const std::string &method, hamiltonian<> &hamiltonian, const std::vector<std::shared_ptr<polariton>> &ppol);
  

} // namespace gmb

#endif //GMB_INIT_H
