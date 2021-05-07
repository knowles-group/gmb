#ifndef GMB_CCSD_EOM_IW_OVOV_H
#define GMB_CCSD_EOM_IW_OVOV_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief W_ovov intermediate
 * 
 */
container<4,double> eom_ccsd_iw_ovov(
    container<2, double> &t1,       ///< CCSD T1 amplitude
    container<4, double> &t2,       ///< tau intermediate
    container<4, double> &i_ooov,   ///> anti-symmetrized integral <ij||ka>
    container<4, double> &i_oovv,   ///> anti-symmetrized integral <ij||ab>
    container<4, double> &i_ovov,   ///> anti-symmetrized integral <ia||jb>
    container<4, double> &i_ovvv    ///> anti-symmetrized integral <ia||bc>

) {

    container<4, double> iw_ovov(i_ovov.get_space());

    libtensor::letter i,j,m,n,a,b,e,f;

    iw_ovov(j|a|i|b) = i_ovov(j|a|i|b)
                     + contract(m, t1(m|a), i_ooov(m|j|i|b))
                     - contract(e, t1(i|e), i_ovvv(j|a|b|e))
                     + contract(m|e, 
                        0.5 * t2(i|m|a|e) + t1(i|e) * t1(m|a), i_oovv(m|j|e|b));
    return iw_ovov;

}

#endif // GMB_CCSD_EOM_IW_OVOV_H
