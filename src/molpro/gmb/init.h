#ifndef GMB_INIT_H
#define GMB_INIT_H

// #include "get_integral.h"
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
  void init(const std::string &filename, const std::string &method, hamiltonian<> &hamiltonian, const std::vector<std::shared_ptr<polariton>> &ppol);


  // read fcidump filee
  void read_dump(const std::string &filename, 
               const std::vector<std::shared_ptr<polariton>> &v_ppol,
               std::vector<std::vector<bool>>& v_exist,
               std::vector<std::vector<size_t>>& v_norb,
               const std::vector<orb_type>& v_orb_type, 
               std::vector<std::vector<std::pair<syms_t, syms_t>>>& v_psi, 
               std::vector<std::vector<std::vector<int>>>& v_shift,
               std::vector<libtensor::bispace<1>> &v_sp,
               const std::vector<spin>& v_spin,
               bool &uhf);


  container<2,double> get_integral(const std::string &fname_integrals, const std::string &fname_header, 
    const std::vector<std::shared_ptr<polariton>> &v_ppol, const orb_type &o1, const orb_type &o2, bool add_ph = true);
  
  container<4,double> get_integral(const std::string &filename, 
    const std::vector<std::shared_ptr<polariton>> &v_ppol,
    const orb_type &o1, const orb_type &o2, const orb_type &o3, const orb_type &o4);
  
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

  container<4,double> get_i(const std::string &filename, 
               const std::vector<std::shared_ptr<polariton>> &v_ppol,
               const orb_type &o1, const orb_type &o2, const orb_type &o3, const orb_type &o4);

  void get_one_electron_part(container<2,double> &integral, 
               const std::string &filename, 
               const std::vector<std::vector<bool>> &v_exist,
               const std::vector<std::vector<size_t>>& v_norb,
               const std::vector<orb_type> &v_orb_type, 
               const std::vector<std::vector<std::pair<syms_t, syms_t>>>& v_psi,
               const std::vector<std::vector<std::vector<int>>>& v_shift, 
               const bool &uhf);

  void get_one_photon_part(container<2,double> &integral, 
               const std::vector<std::shared_ptr<polariton>> &v_ppol,
               const std::vector<std::vector<bool>>& v_exist,
               const std::vector<orb_type>& v_orb_type);
  
  void get_two_electron_part(container<2,double> &integral, 
               const std::string &filename, 
               const std::vector<std::vector<bool>> &v_exist,
               const std::vector<std::vector<size_t>>& v_norb,
               const std::vector<orb_type> &v_orb_type, 
               const std::vector<std::vector<std::pair<syms_t, syms_t>>>& v_psi,
               const std::vector<std::vector<std::vector<int>>>& v_shift, 
               const bool &uhf);

 void get_electron_photon_part(container<4,double> &integral, 
               const std::vector<std::shared_ptr<polariton>> &v_ppol,
               const std::vector<std::vector<bool>> &v_exist,
               const std::vector<std::vector<size_t>>& v_norb,
               const std::vector<orb_type> &v_orb_type, 
               const std::vector<std::vector<std::pair<syms_t, syms_t>>>& v_psi,
               const std::vector<std::vector<std::vector<int>>>& v_shift);

 
  container<2,double> set_space(const std::vector<orb_type> &v_orb_type, const std::vector<libtensor::bispace<1>> &v_sp);
              
} // namespace gmb

#endif //GMB_INIT_H
