#ifndef GMB_GET_INTEGRAL_H
#define GMB_GET_INTEGRAL_H


#include "molpro/FCIdump.h"
#include "libtensor/libtensor.h"
#include "libtensor_utils.h"
#include "container.h"
#include <vector>
#include <string>
#include "get_integral.h"

using sym_t=unsigned int;
using syms_t=std::vector<sym_t>; ///< position in orbital space for each symmetry

// orbital types occupied, virtual & basis (
enum orb_type{o, v, b};
enum spin{alpha=0, beta=1, photon=2};

struct polariton {
  sym_t nmax{0};
  double gamma{0.0};
  double omega{0.0};
  double lambd{0.0};
  polariton() = default;
  polariton(sym_t nmax_, double gamma_, double omega_) 
  : nmax{nmax_}, gamma{gamma_}, omega{omega_} 
  {lambd = gamma*(sqrt(2*omega));}
};

// get nuclear energy
double get_integral(std::string filename);

// read fcidump file
void read_dump(std::string filename, 
               const std::vector<orb_type>& orb_types, 
               const std::vector<spin>& v_spin,
               std::vector<std::vector<std::pair<syms_t, syms_t>>>& v_psi, 
               std::vector<std::vector<size_t>>& v_norb,
               std::vector<std::vector<std::vector<int>>>& v_shift,
               std::vector<libtensor::bispace<1>>& v_sp,
               std::vector<std::vector<bool>>& ssss);
               
// get one-electron integral
container<2,double> get_integral(std::string filename, 
                                 orb_type o1, 
                                 orb_type o2,
                                 bool so_basis = true);

// get two-electron integral
container<4,double> get_integral(std::string filename, 
                                 orb_type o1, 
                                 orb_type o2, 
                                 orb_type o3, 
                                 orb_type o4,
                                 bool so_basis = true);

// get anti-symmetrized two-electron integral
container<4,double> get_i(std::string filename, 
                          orb_type o1, 
                          orb_type o2, 
                          orb_type o3, 
                          orb_type o4);
#endif //GMB_GET_INTEGRAL_H