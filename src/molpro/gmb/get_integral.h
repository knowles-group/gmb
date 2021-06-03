#ifndef GMB_GET_INTEGRAL_H
#define GMB_GET_INTEGRAL_H


#include "molpro/FCIdump.h"
#include "libtensor/libtensor.h"
#include "libtensor_utils.h"
#include "container.h"
#include <vector>
#include <string>


using sym_t=unsigned int;
using syms_t=std::vector<sym_t>; ///< position in orbital space for each symmetry

// orbital types occupied, virtual & basis (
enum orb_type{o, v, b};
enum spin{alpha=0, beta=1, photon=2};

// get nuclear energy
double get_integral(std::string filename);

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