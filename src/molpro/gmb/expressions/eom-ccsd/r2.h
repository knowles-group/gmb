#ifndef GMB_CCSD_R2_H
#define GMB_CCSD_R2_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief Action to be applied to EOMEE-CCSD R2 amplitudes.
 * 
 */
container<4,double> eom_ccsd_r2 (
    container<2, double> &r1,      ///< EOM-CCSD R1 amplitude
    container<4, double> &r2,      ///< EOM-CCSD R2 amplitude
    container<4, double> &t2,      ///< CCSD T2 amplitude
    container<2, double> &if_oo,   ///< intermediate - oo block
    container<2, double> &if_vv,   ///< intermediate - vv block
    container<4, double> &i_oovv,  ///> anti-symmetrized integral <ij||ab>
    container<4, double> &iw_oooo, ///> intermediate - oooo block
    container<4, double> &iw2_ooov, ///> intermediate 2 - ooov block
    container<4, double> &iw_ovov, ///> intermediate - ovov block
    container<4, double> &iw_ovvv,  ///> intermediate - ovvv block
    container<4, double> &iw2_ovvv, ///> intermediate - ovvv block
    container<4, double> &iw_vvvv  ///> intermediate - vvvv block
) {

    container<4, double> r2_new(t2.get_space());

    libtensor::letter i,j,m,n,a,b,e,f;
        r2_new(i|j|a|b) = asymm(a, b, contract(e, r2(i|j|a|e), if_vv(b|e)))
                   - asymm(i, j, contract(m, r2(i|m|a|b), if_oo(m|j)))
                   - asymm(a, b, contract(m, r1(m|a), iw2_ooov(i|j|m|b)))
                   + asymm(i, j, contract(e, r1(i|e), iw_ovvv(j|e|b|a)))
                   + 0.5 * contract(m|n, r2(m|n|a|b), iw_oooo(m|n|i|j))
                   + asymm(a, b, asymm(i, j, 
                        contract(e|m, r2(i|m|a|e), iw_ovov(m|b|j|e))))
                   + 0.5 * contract(e|f, r2(i|j|e|f), iw_vvvv(a|b|e|f))
                   + asymm(a, b, contract(m|f, r1(m|f), 
                        contract(e, t2(i|j|a|e), iw2_ovvv(m|b|f|e))))
                   - asymm(i, j, contract(n|e, r1(n|e), 
                        contract(m, t2(i|m|a|b), iw2_ooov(m|n|j|e))))
                   - asymm(i, j, contract(n|e|f, r2(j|n|e|f), 
                        contract(m, t2(i|m|a|b), i_oovv(m|n|e|f))))
                   - asymm(a, b, contract(m|n|f, r2(m|n|b|f), 
                        contract(e, t2(i|j|a|e), i_oovv(m|n|e|f))));

    return r2_new;
}     

#endif // GMB_CCSD_R2_H
