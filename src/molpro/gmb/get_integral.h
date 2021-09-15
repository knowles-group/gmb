#ifndef GMB_GET_INTEGRAL_H
#define GMB_GET_INTEGRAL_H


#include "molpro/FCIdump.h"
#include <libtensor/libtensor.h>
#include "libtensor_utils.h"
#include "container.h"
#include <vector>
#include <string>

using sym_t=unsigned int;
using syms_t=std::vector<sym_t>; ///< position in orbital space for each symmetry

// orbital types occupied, virtual & basis (
enum orb_type{o, v, b};
// possible spins (
enum spin{alpha=0, beta=1, photon=2};

struct polariton {
  sym_t nmax{0};
  double gamma{0.0};
  double omega{0.0};
  double lambd{0.0};
  std::string fname_dm;
  std::string fname_sm;

  polariton() = default;
  polariton(sym_t nmax_, double gamma_, double omega_) 
  : nmax{nmax_}, gamma{gamma_}, omega{omega_} 
  {
    lambd = gamma*(sqrt(2*omega));
  }

};

// get nuclear energy
double get_integral(const std::string &filename);

/**
 * @brief Read fcidump file
 * 
 * @param filename 
 * @param v_ppol pointers to polariton modes
 * @param v_exist if block exists or not
 * @param v_norb number of orbitals in each bra/ket
 * @param v_orb_type orbital types
 * @param v_psi bra and ket
 * @param v_shift symmetry shift
 * @param v_sp 1D spaces for each bra/ket
 * @param v_spin possible spins
 * @param uhf 
 */
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
               const orb_type &o1, const orb_type &o2, const orb_type &o3, const orb_type &o4, const bool &add_ph = true);

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
              


#endif //GMB_GET_INTEGRAL_H