#ifndef GMB_CCSD_T2_H
#define GMB_CCSD_T2_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief Action to be applied to CCSD T2 amplitudes.
 * 
 */
container<4,double> ccsd_t2 (
    container<2, double> &t1,      ///< CCSD T1 amplitude
    container<4, double> &t2,      ///< CCSD T2 amplitude
    container<4, double> &tau,     ///< tau intermediate
    container<2, double> &if_oo,   ///< intermediate - oo block
    container<2, double> &if_vv,   ///< intermediate - vv block
    container<4, double> &i_ooov,  ///> anti-symmetrized integral <ij||ka>
    container<4, double> &i_oovv,  ///> anti-symmetrized integral <ij||ab>
    container<4, double> &i_ovov,  ///> anti-symmetrized integral <ia||jb>
    container<4, double> &i_ovvv,  ///> anti-symmetrized integral <ia||bc>
    container<4, double> &iw_oooo, ///> intermediate - oooo block
    container<4, double> &iw_ovov, ///> intermediate - ovov block
    container<4, double> &iw_vvvv  ///> intermediate - vvvv block
) {

    container<4, double> t2_new(t2.get_space());

    libtensor::letter i,j,m,n,a,b,e,f;
        t2_new(i|j|a|b) = i_oovv(i|j|a|b)
                   + asymm(a, b, contract(e, t2(i|j|a|e), if_vv(b|e)))
                   - asymm(i, j, contract(m, t2(i|m|a|b), if_oo(j|m)))
                   + 0.5 * contract(m|n, tau(m|n|a|b), iw_oooo(m|n|i|j))
                   + 0.5 * contract(e|f, tau(i|j|e|f), iw_vvvv(a|b|e|f))
                   + asymm(a, b, asymm(i, j, 
                        contract(e|m, t2(i|m|a|e), iw_ovov(m|b|j|e))
                      - contract(e, t1(i|e), 
                            contract(m, t1(m|b), i_ovov(j|e|m|a)))))
                   + asymm(i, j, contract(e, t1(i|e), i_ovvv(j|e|b|a)))
                   - asymm(a, b, contract(m, t1(m|a), i_ooov(i|j|m|b)));

    return t2_new;
}     

#endif // GMB_CCSD_T2_H
