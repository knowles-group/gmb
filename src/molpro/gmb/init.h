#ifndef GMB_INIT_H
#define GMB_INIT_H

#include "get_integral.h"
#include "utils.h"

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
  void init(const std::string &filename, const std::string &method, hamiltonian<> &hamiltonian, const std::vector<std::unique_ptr<polariton>> &ppol);


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

  // read fcidump filee
  void read_dump(const std::string &filename, 
                 const std::vector<std::unique_ptr<polariton>> &v_ppol,
                 const std::vector<orb_type>& orb_types, 
                 const std::vector<spin>& v_spin,
                 std::vector<std::vector<std::pair<syms_t, syms_t>>>& v_psi, 
                 std::vector<std::vector<size_t>>& v_norb,
                 std::vector<std::vector<std::vector<int>>>& v_shift,
                 std::vector<libtensor::bispace<1>>& v_space,
                 std::vector<std::vector<bool>>& sssv_exist);

  container<2,double> get_hamiltonian(const std::string &filename, 
    const std::vector<std::unique_ptr<polariton>> &v_ppol, const orb_type &o1, const orb_type &o2);


} // namespace gmb

#endif //GMB_INIT_H
