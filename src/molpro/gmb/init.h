#ifndef GMB_SRC_MOLPRO_GMB_INIT_H
#define GMB_SRC_MOLPRO_GMB_INIT_H

#include "get_integral.h"
#include "utils.h"

#include "hamiltonian.h"
#include "supercontainer.h"

#include <memory>
#include <numeric>

#include <molpro/iostream.h>

namespace gmb {

  /**
   * @brief Initialise Hamiltonian 
   * 
   * @param filename 
   * @param method 
   * @param hamiltonian 
   */
  void init(const std::string &filename, 
            const std::string &method, 
            hamiltonian<> &hamiltonian, 
            const std::vector<std::shared_ptr<polariton>> &ppol, 
            const std::vector<std::unique_ptr<vibration>> &pvib);

  /**
   * @brief Read fock matrix for orbital energies
   * 
   * @param f_xx block of fock matrix to read 
   * @param ss where to print
   * @param x O or V for Occupied/Virtual
   */
  void readf(container<2> &f_xx, std::ostringstream &ss, const char &x, 
            const std::vector<std::shared_ptr<polariton>> &v_ppol); 

} // namespace gmb

#endif // GMB_SRC_MOLPRO_GMB_INIT_H
