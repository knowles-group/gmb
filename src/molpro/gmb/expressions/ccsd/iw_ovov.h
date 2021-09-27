#ifndef GMB_SRC_MOLPRO_GMB_CCSD_IW_OVOV_H
#define GMB_SRC_MOLPRO_GMB_CCSD_IW_OVOV_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief W_ovov intermediate
 * 
 */
container<4,double> ccsd_iw_ovov(
    container<2, double> &t1,       ///< CCSD T1 amplitude
    container<4, double> &t2,       ///< tau intermediate
    container<4, double> &i_ooov,   ///> anti-symmetrized integral <ij||ka>
    container<4, double> &i_oovv,   ///> anti-symmetrized integral <ij||ab>
    container<4, double> &i_ovov,   ///> anti-symmetrized integral <ia||jb>
    container<4, double> &i_ovvv    ///> anti-symmetrized integral <ia||bc>

) {

    container<4, double> iw_ovov(i_ovov.get_space());

    libtensor::letter i,j,m,n,a,b,e,f;

    iw_ovov(m|b|j|e) = - i_ovov(m|b|j|e)
                     - contract(n, t1(n|b), i_ooov(n|m|j|e))
                     + contract(f, t1(j|f), i_ovvv(m|b|e|f))
                     + contract(n|f, 
                        0.5 * t2(j|n|b|f) - t1(j|f) * t1(n|b), i_oovv(m|n|e|f));
    return iw_ovov;

}

#endif // GMB_SRC_MOLPRO_GMB_CCSD_IW_OVOV_H
