#ifndef GMB_INIT_H
#define GMB_INIT_H

#include "get_integral.h"

#include "hamiltonian.h"
#include "supercontainer.h"

namespace gmb {

  /**
   * @brief Initialise Hamiltonian 
   * 
   * @param filename 
   * @param method 
   * @param hamiltonian 
   */
  void init(const std::string &filename, const std::string &method, hamiltonian<> &hamiltonian, const std::unique_ptr<polariton> &ppol);


  /**
   * @brief get anti-symmetrized two-electron integral
   * 
   * @param filename 
   * @param o1 
   * @param o2 
   * @param o3 
   * @param o4 
   * @return container<4,double> 
   */
  container<4,double> get_i(std::string filename, 
                            orb_type o1, 
                            orb_type o2, 
                            orb_type o3, 
                            orb_type o4);


  void add_self_energy(std::string method, hamiltonian<> &hamiltonian);

} // namespace gmb

#endif //GMB_INIT_H
