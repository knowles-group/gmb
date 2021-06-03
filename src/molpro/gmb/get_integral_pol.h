#ifndef GMB_GET_INTEGRAL_POL_H
#define GMB_GET_INTEGRAL_POL_H


#include "molpro/FCIdump.h"
#include "libtensor/libtensor.h"
#include "libtensor_utils.h"
#include "container.h"
#include <vector>
#include <string>
#include "get_integral.h"


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