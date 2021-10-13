#ifndef GMB_SRC_MOLPRO_GMB_EOM_CCSD_L1_H
#define GMB_SRC_MOLPRO_GMB_EOM_CCSD_L1_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief Action to be applied to EOMEE-CCSD L1 amplitudes.
 * 
 */
container<2,double> eom_ccsd_l1(
    container<2, double> &l1,         ///< EOM-CCSD L1 amplitude
    container<4, double> &l2,         ///< EOM-CCSD L2 amplitude
    container<4, double> &t2,         ///< CCSD T2 amplitude
    container<2, double> &il_oo,      ///< intermediate
    container<2, double> &il_vv,      ///< intermediate
    container<2, double> &if_oo,      ///< intermediate - oo block
    container<2, double> &if_ov,      ///< intermediate - ov block
    container<2, double> &if_vv,      ///< intermediate - vv block
    container<4, double> &iw_ovov,    ///> anti-symmetrized integral <ij||ka>
    container<4, double> &iw_ooov,   ///> anti-symmetrized integral <ia||jb>
    container<4, double> &iw2_ooov,   ///> anti-symmetrized integral <ia||jb>
    container<4, double> &iw_ovvv,   ///> anti-symmetrized integral <ia||bc>
    container<4, double> &iw2_ovvv,   ///> anti-symmetrized integral <ia||bc>
    double l0 = 0.0                        ///< EOM-CCSD L0 

) {

    container<2, double> l1_new(l1.get_space());

    libtensor::letter i,j,m,n,o,a,b,e,f,g;

    l1_new(i|a) = 
                l0*if_ov(i|a)
                + contract(e, l1(i|e), if_vv(e|a))
                - contract(m, l1(m|a), if_oo(i|m))
                - contract(e|m, l1(m|e), iw_ovov(m|a|i|e))  
                - 0.5 * contract(e|m|n, l2(m|n|a|e), iw_ooov(m|n|i|e))
                + 0.5 * contract(e|f|m, l2(i|m|e|f), iw_ovvv(m|a|f|e))
                + contract(e|f, iw2_ovvv(i|e|a|f), il_vv(e|f))
                + contract(m|n, iw2_ooov(i|n|m|a), il_oo(m|n));

    return l1_new;

}

#endif // GMB_SRC_MOLPRO_GMB_EOM_CCSD_L1_H
