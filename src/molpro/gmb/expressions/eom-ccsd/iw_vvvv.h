#ifndef GMB_EOM_CCSD_IW_VVVV_H
#define GMB_EOM_CCSD_IW_VVVV_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief W_vvvv intermediate
 * 
 */
container<4,double> eom_ccsd_iw_vvvv(
    container<2, double> &t1,       ///< CCSD T1
    container<4, double> &tau,      ///< tau intermediate
    container<4, double> &i_oovv,   ///> anti-symmetrized integral <ij||ab>
    container<4, double> &i_ovvv,   ///> anti-symmetrized integral <ia||bc>
    container<4, double> &i_vvvv    ///> anti-symmetrized integral <ab||cd>

) {

    container<4, double> iw_vvvv(i_vvvv.get_space());

    libtensor::letter m,n,a,b,c,d;

    iw_vvvv(c|d|a|b) = i_vvvv(c|d|a|b)
                     + asymm(c, d, 
                       contract(m, t1(m|d), i_ovvv(m|c|a|b)))
                     + 0.5 * contract(m|n, tau(m|n|c|d), i_oovv(m|n|a|b));
                  
    return iw_vvvv;

}


#endif // GMB_EOM_CCSD_IW_VVVV_H
