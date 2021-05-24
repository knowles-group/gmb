#ifndef GMB_EOM_CCSD_IF_VV_H
#define GMB_EOM_CCSD_IF_VV_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief F_vv intermediate
 * 
 */
container<2,double> eom_ccsd_if_vv(
    container<2, double> &t1,         ///< CCSD T1
    container<4, double> &t2,         ///< CCSD T2
    container<2, double> &f_ov,       ///< fock matrix - ov block
    container<2, double> &f_vv,       ///< fock matrix - vv block
    container<4, double> &i_oovv,     ///> anti-symmetrized integral <ij||ab>
    container<4, double> &i_ovvv      ///> anti-symmetrized integral <ia||bc>
) {

    container<2, double> if_vv(f_vv.get_space());

    libtensor::letter m,n,a,b,e;

    if_vv(a|b) = f_vv(a|b)
               - contract(m, t1(m|a), f_ov(m|b))
               + contract(e|m, t1(m|e), i_ovvv(m|a|e|b))
               - contract(e|m|n, t1(m|a)*t1(n|e), i_oovv(m|n|b|e))
               - 0.5 * contract(e|m|n, t2(m|n|a|e), i_oovv(m|n|b|e));
    return if_vv;

}

#endif //GMB_EOM_CCSD_IF_VV_H
