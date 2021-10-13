#ifndef GMB_SRC_MOLPRO_GMB_CCSD_L2_H
#define GMB_SRC_MOLPRO_GMB_CCSD_L2_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief Action to be applied to EOMEE-CCSD R2s.
 * 
 */
container<4,double> eom_ccsd_l2 (
    container<2, double> &l1,       ///< EOM-CCSD R1
    container<4, double> &l2,       ///< EOM-CCSD R2
    container<4, double> &t2,       ///< CCSD T2
    container<2, double> &if_oo,    ///< intermediate if_oo 
    container<2, double> &if_ov,    ///< intermediate if_ov 
    container<2, double> &if_vv,    ///< intermediate if_vv 
    container<2, double> &il_oo,   ///< intermediate il_oo
    container<2, double> &il_vv,   ///< intermediate il_vv
    container<4, double> &i_oovv,   ///> anti-symmetrized integral <ij||ab>
    container<4, double> &iw_oooo,  ///> intermediate iw_oooo 
    container<4, double> &iw_ooov,  ///> intermediate iw_ooov 
    container<4, double> &iw2_ooov, ///> intermediate iw2_ooov 
    container<4, double> &iw_ovov,  ///> intermediate iw_ovov 
    container<4, double> &iw_ovvv,  ///> intermediate iw_ovvv 
    container<4, double> &iw2_ovvv, ///> intermediate iw2_ovvv 
    container<4, double> &iw_vvvv,   ///> intermediate iw_vvvv 
    double l0 = 0.0                     ///< EOM-CCSD L0 
) {

    container<4, double> l2_new(t2.get_space());

    libtensor::letter i,j,k,l,m,n,a,b,c,d,e,f;
    
    l2_new(i|j|a|b) = 
                    l0*i_oovv(i|j|a|b)
                      + asymm(a, b, 
                          contract(e, l2(i|j|a|e), if_vv(e|b))
                        - contract(m, l1(m|a), iw2_ooov(i|j|m|b))
                        - 0.5*contract(e, i_oovv(i|j|e|b), 
                            contract(m|n|f, t2(m|n|e|f), l2(m|n|a|f))))
                     + asymm(i, j, 
                           - contract(m, l2(i|m|a|b), if_oo(j|m))
                           + contract(e, l1(i|e), iw2_ovvv(j|e|b|a))
                        - 0.5*contract(m, i_oovv(m|j|a|b), 
                            contract(n|e|f, t2(m|n|e|f), l2(i|n|e|f)))
                        + asymm(a,b,
                            l1(i|a)*if_ov(j|b)
                         - contract(e|m, l2(i|m|a|e), iw_ovov(m|b|j|e))))
                     + 0.5 * contract(m|n, l2(m|n|a|b), iw_oooo(i|j|m|n))
                     + 0.5 * contract(e|f, l2(i|j|e|f), iw_vvvv(e|f|a|b));

    return l2_new;
}     

#endif // GMB_SRC_MOLPRO_GMB_CCSD_L2_H
