#ifndef GMB_EOM_CCSD_IW_OOOO_H
#define GMB_EOM_CCSD_IW_OOOO_H

#include <libtensor/libtensor.h>
#include "../../container.h"


/**
 * @brief W_oooo intermediate
 * 
 */
container<4,double> eom_ccsd_iw_oooo(
    container<2, double> &t1,      ///< CCSD T1
    container<4, double> &tau,     ///< tau intermediate
    container<4, double> &i_oooo,  ///> anti-symmetrized integral <ij||kl>
    container<4, double> &i_ooov,  ///> anti-symmetrized integral <ij||ka>
    container<4, double> &i_oovv   ///> anti-symmetrized integral <ij||ab>

) {

    container<4, double> iw_oooo(i_oooo.get_space());

    libtensor::letter i,j,k,l,e,f;

    iw_oooo(i|j|k|l) = i_oooo(i|j|k|l)
                     + asymm(i, j, 
                        contract(e, t1(j|e), i_ooov(k|l|i|e)))
                     + 0.5 * contract(e|f, tau(i|j|e|f), i_oovv(k|l|e|f));

    return iw_oooo;

}


#endif //GMB_EOM_CCSD_IW_OOOO_H
