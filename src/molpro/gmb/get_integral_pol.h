#ifndef GMB_GET_INTEGRAL_POL_H
#define GMB_GET_INTEGRAL_POL_H


#include "molpro/FCIdump.h"
#include "libtensor/libtensor.h"
#include "libtensor_utils.h"
#include "container.h"
#include <vector>
#include <string>
#include "get_integral.h"

void read_dump(std::string filename, 
               const std::vector<orb_type>& orb_types, 
               const std::vector<spin>& v_spin,
               std::vector<std::vector<std::pair<syms_t, syms_t>>>& v_psi, 
               std::vector<std::vector<size_t>>& v_norb,
               std::vector<std::vector<std::vector<int>>>& v_shift,
               std::vector<libtensor::bispace<1>>& v_sp,
               std::vector<std::vector<bool>>& ssss);
               
// get one-electron integral
container<2,double> get_integral_pol(std::string filename, 
                                 orb_type o1, 
                                 orb_type o2);

// get two-electron integral
container<4,double> get_integral_pol(std::string filename, 
                                 orb_type o1, 
                                 orb_type o2, 
                                 orb_type o3, 
                                 orb_type o4,
                                 bool add_ph = false);

// get anti-symmetrized two-electron integral
container<4,double> get_i_pol(std::string filename, 
                          orb_type o1, 
                          orb_type o2, 
                          orb_type o3, 
                          orb_type o4);
#endif //GMB_GET_INTEGRAL_POL_H