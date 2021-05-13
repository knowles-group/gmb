#ifndef GMB_EOM_CCSD_IR1_VV_H
#define GMB_EOM_CCSD_IR1_VV_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief W_ovvv intermediate
 * 
 */
container<2,double> eom_ccsd_ir1_vv(
    container<2, double> &f_vv,     ///< fock operator - vv block
    container<2, double> &r1,       ///< EOM-CCSD R1
    container<4, double> &iw2_ovvv  ///> intermediate 

) {

    container<2, double> ir1_vv(f_vv.get_space());

    libtensor::letter i,j,m,n,a,b,c,e,f;

    ir1_vv(a|e) = contract(m|f, r1(m|f), iw2_ovvv(m|a|e|f))
    return ir1_vv;

}

#endif //GMB_EOM_CCSD_IR1_VV_H
