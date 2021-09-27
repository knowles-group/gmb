#ifndef GMB_SRC_MOLPRO_GMB_EOM_CCSD_IF_OV_H
#define GMB_SRC_MOLPRO_GMB_EOM_CCSD_IF_OV_H

#include <libtensor/libtensor.h>
#include "../../container.h"


/**
 * @brief F_ov intermediate
 * 
 */
container<2,double> eom_ccsd_if_ov(
    container<2, double> &t1,     ///< CCSD T1
    container<2, double> &f_ov,   ///< fock matrix - ov block
    container<4, double> &i_oovv  ///> anti-symmetrized integral <ij||ab>

) {

    container<2, double> if_ov(f_ov.get_space());

    libtensor::letter i,m,a,e;

    if_ov(i|a) = f_ov(i|a)
               + contract(m|e, t1(m|e), i_oovv(i|m|a|e));
    return if_ov;

}


#endif // GMB_SRC_MOLPRO_GMB_EOM_CCSD_IF_OV_H
