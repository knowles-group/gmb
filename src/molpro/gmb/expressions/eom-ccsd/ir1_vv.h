#ifndef GMB_SRC_MOLPRO_GMB_EOM_CCSD_IR1_VV_H
#define GMB_SRC_MOLPRO_GMB_EOM_CCSD_IR1_VV_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief ir1_vv intermediate
 * 
 */
container<2,double> eom_ccsd_ir1_vv(
    container<2, double> &f_vv,     ///< fock operator - vv block
    container<2, double> &r1,       ///< EOM-CCSD R1 
    container<4, double> &iw2_ovvv  ///> intermediate 

) {

    container<2, double> ir1_vv(f_vv.get_space());

    libtensor::letter m,b,e,f;

    ir1_vv(b|e) = contract(m|f, r1(m|f), iw2_ovvv(m|b|f|e));
    return ir1_vv;

}

#endif // GMB_SRC_MOLPRO_GMB_EOM_CCSD_IR1_VV_H
