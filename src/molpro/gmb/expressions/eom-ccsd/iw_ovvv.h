#ifndef GMB_EOM_CCSD_IW_OVVV_H
#define GMB_EOM_CCSD_IW_OVVV_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief W_ovvv intermediate
 * 
 */
container<4,double> eom_ccsd_iw_ovvv(
    container<2, double> &t1,       ///< CCSD T1 amplitude
    container<4, double> &t2,       ///< CCSD T2 amplitude
    container<4, double> &tau,      ///< tau intermediate
    container<2, double> &f_ov,     ///< fock matrix - ov block
    container<4, double> &i_ooov,   ///> anti-symmetrized integral <ij||ka>
    container<4, double> &i_oovv,   ///> anti-symmetrized integral <ij||ab>
    container<4, double> &i_ovov,   ///> anti-symmetrized integral <ia||kb>
    container<4, double> &i_ovvv,   ///> anti-symmetrized integral <ia||jb>
    container<4, double> &i_vvvv    ///> anti-symmetrized integral <ab||cd>

) {

    container<4, double> iw_ovvv(i_ovvv.get_space());

    libtensor::letter i,j,m,n,a,b,c,e,f;

    iw_ovvv(i|a|b|c) = i_ovvv(i|a|b|c)
                     - contract(m, t2(i|m|c|b), f_ov(m|a))
                     + contract(e, t1(i|e), i_vvvv(b|c|a|e))
                     + asymm(b, c,
                        - contract(m, t1(m|c), i_ovov(m|b|i|a))
                        + contract(m|e, t2(i|m|c|e), i_ovvv(m|b|e|a))
                        - contract(m|e, t1(i|e)*t1(m|c), i_ovvv(m|b|a|e))
                        - contract(m, t1(m|b), 
                            contract(e|n, t2(n|i|e|c), i_oovv(m|n|a|e))))
                     + 0.5 * contract(m|n, tau(m|n|c|b), i_ooov(m|n|i|a))
                     - contract(m|e, t1(m|e), 
                        contract(n, t2(n|i|b|c), i_oovv(n|m|a|e)))
                     + 0.5 * contract(m|n, tau(m|n|c|b), 
                        contract(e, t1(i|e), i_oovv(m|n|e|a)));
    return iw_ovvv;

}

#endif // GMB_EOM_CCSD_IW_OVVV_H
