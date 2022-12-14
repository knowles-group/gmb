#ifndef GMB_SRC_MOLPRO_GMB_CCSD_R2_H
#define GMB_SRC_MOLPRO_GMB_CCSD_R2_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief Action to be applied to EOMEE-CCSD R2s.
 * 
 */
container<4,double> eom_ccsd_r2 (
    container<2, double> &r1,       ///< EOM-CCSD R1
    container<4, double> &r2,       ///< EOM-CCSD R2
    container<4, double> &t2,       ///< CCSD T2
    container<2, double> &if_oo,    ///< intermediate if_oo 
    container<2, double> &if_vv,    ///< intermediate if_vv 
    container<2, double> &ir1_oo,   ///< intermediate ir1_oo
    container<2, double> &ir2_oo,   ///< intermediate ir2_oo
    container<2, double> &ir1_vv,   ///< intermediate ir1_vv
    container<2, double> &ir2_vv,   ///< intermediate ir2_vv
    container<4, double> &i_oovv,   ///> anti-symmetrized integral <ij||ab>
    container<4, double> &iw_oooo,  ///> intermediate iw_oooo 
    container<4, double> &iw_ooov,  ///> intermediate iw_ooov 
    container<4, double> &iw2_ooov, ///> intermediate iw2_ooov 
    container<4, double> &iw_ovov,  ///> intermediate iw_ovov 
    container<4, double> &iw_ovvv,  ///> intermediate iw_ovvv 
    container<4, double> &iw2_ovvv, ///> intermediate iw2_ovvv 
    container<4, double> &iw_vvvv   ///> intermediate iw_vvvv 
) {

    container<4, double> r2_new(t2.get_space());

    libtensor::letter i,j,k,l,m,n,a,b,c,d,e,f;
    r2_new(i|j|a|b) = asymm(a, b, 
                          contract(e, r2(i|j|a|e), if_vv(b|e))
                        - contract(m, r1(m|a), iw_ooov(i|j|m|b))
                        + contract(e, t2(i|j|a|e), ir1_vv(b|e))
                        - contract(e, t2(i|j|a|e), ir2_vv(b|e)))
                     + asymm(i, j, 
                           - contract(m, r2(i|m|a|b), if_oo(m|j))
                           + contract(e, r1(i|e), iw_ovvv(j|e|b|a))
                           - contract(m, t2(i|m|a|b), ir1_oo(j|m))
                           - contract(m, t2(i|m|a|b), ir2_oo(j|m))
                        + asymm(a, b, 
                         - contract(e|m, r2(i|m|a|e), iw_ovov(m|b|j|e))))
                     + 0.5 * contract(m|n, r2(m|n|a|b), iw_oooo(i|j|m|n))
                     + 0.5 * contract(e|f, r2(i|j|e|f), iw_vvvv(a|b|e|f));
    return r2_new;
}     

#endif // GMB_SRC_MOLPRO_GMB_CCSD_R2_H
