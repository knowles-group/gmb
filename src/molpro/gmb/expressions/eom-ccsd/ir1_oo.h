#ifndef GMB_EOM_CCSD_IR1_OO_H
#define GMB_EOM_CCSD_IR1_OO_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief ir1_oo intermediate
 * 
 */
container<2,double> eom_ccsd_ir1_oo(
    container<2, double> &f_oo,     ///< fock operator - vv block
    container<2, double> &r1,       ///< EOM-CCSD R1 
    container<4, double> &iw2_ooov  ///> intermediate iw2_ooov

) {

    container<2, double> ir1_oo(f_oo.get_space());

    libtensor::letter i,j,m,n,a,b,c,e,f;

    ir1_oo(j|m) = contract(n|e, r1(n|e), iw2_ooov(m|n|j|e));
    return ir1_oo;

}

#endif //GMB_EOM_CCSD_IR1_OO_H
