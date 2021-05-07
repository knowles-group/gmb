#ifndef GMB_EOM_CCSD_IW_OOOV_H
#define GMB_EOM_CCSD_IW_OOOV_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief W_ooov intermediate
 * 
 */
container<4,double> eom_ccsd_iw_ooov(
    container<2, double> &t1,       ///< CCSD T1 amplitude
    container<4, double> &t2,       ///< CCSD T2 intermediate
    container<4, double> &tau,      ///< tau intermediate
    container<2, double> &f_ov,     ///< fock matrix - ov block
    container<4, double> &i_oooo,   ///> anti-symmetrized integral <ij||kl>
    container<4, double> &i_ooov,   ///> anti-symmetrized integral <ij||ka>
    container<4, double> &i_oovv,   ///> anti-symmetrized integral <ij||ab>
    container<4, double> &i_ovov,   ///> anti-symmetrized integral <ia||jb>
    container<4, double> &i_ovvv    ///> anti-symmetrized integral <ia||jb>

) {

    container<4, double> iw_ooov(i_ooov.get_space());

    libtensor::letter i,j,k,m,n,a,b,e,f;

    iw_ooov(i|j|k|a) = i_ooov(i|j|k|a)
                     + contract(e, t2(j|i|a|e), f_ov(k|e))
                     + asymm(i, j,
                        contract(e, t1(j|e), i_ovov(k|a|i|e)))
                     - contract(m, t1(m|a), i_oooo(m|k|j|i))

                     + asymm(i, j,
                       - contract(m|e, t2(j|m|a|e), i_ooov(m|k|i|e))
                       + contract(m|e, t1(j|e)*t1(m|a), i_ooov(m|k|i|e)))
                     - 0.5 * contract(e|f, tau(j|i|e|f), i_ovvv(k|a|e|f))

                     + contract(m|f, t1(m|f), 
                        contract(e, t2(j|i|a|e), i_oovv(k|m|e|f)))
                     + asymm(i, j,
                        contract(e, t1(i|e), 
                          contract(m|f, t2(j|m|a|f), i_oovv(k|m|e|f))))

                     - 0.5 * contract(m, t1(m|a), 
                           contract(e|f, tau(j|i|e|f), i_oovv(m|k|e|f)))
                     ;


    return iw_ooov;

}

#endif // GMB_EOM_CCSD_IW_OOOV_H
