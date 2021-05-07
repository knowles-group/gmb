#ifndef GMB_EOM_CCSD_R1_H
#define GMB_EOM_CCSD_R1_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief Action to be applied to EOMEE-CCSD R1 amplitudes.
 * 
 */
container<2,double> eom_ccsd_r1(
    container<2, double> &r1,      ///< EOM-CCSD R1 amplitude
    container<4, double> &r2,      ///< EOM-CCSD R2 amplitude
    container<2, double> &if_oo,   ///< intermediate - oo block
    container<2, double> &if_ov,   ///< intermediate - ov block
    container<2, double> &if_vv,   ///< intermediate - vv block
    container<4, double> &iw_ovov,  ///> anti-symmetrized integral <ij||ka>
    container<4, double> &iw2_ooov,  ///> anti-symmetrized integral <ia||jb>
    container<4, double> &iw2_ovvv   ///> anti-symmetrized integral <ia||bc>

) {

    container<2, double> r1_new(r1.get_space());

    libtensor::letter i,j,m,n,a,b,e,f;

    r1_new(i|a) = contract(e, r1(i|e), if_vv(a|e))
                - contract(m, r1(m|a), if_oo(m|i))
                + contract(e|m, r2(i|m|a|e), if_ov(m|e))
                - contract(e|m, r1(m|e), iw_ovov(m|a|i|e))  
                - 0.5 * contract(e|f|m, r2(i|m|e|f), iw2_ovvv(m|a|e|f))
                - 0.5 * contract(e|m|n, r2(m|n|a|e), iw2_ooov(m|n|i|e))
                ;

    return r1_new;

}

#endif // GMB_EOM_CCSD_R1_H
