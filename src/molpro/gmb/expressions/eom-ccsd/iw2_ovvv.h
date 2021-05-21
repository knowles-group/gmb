#ifndef GMB_EOM_CCSD_IW2_OVVV_H
#define GMB_EOM_CCSD_IW2_OVVV_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief W_ovvv intermediate
 * 
 */
container<4,double> eom_ccsd_iw2_ovvv(
    container<2, double> &t1,       ///< CCSD T1
    container<4, double> &i_oovv,   ///> anti-symmetrized integral <ij||ab>
    container<4, double> &i_ovvv    ///> anti-symmetrized integral <ia||bc>

) {

    container<4, double> iw_ovvv(i_ovvv.get_space());

    // libtensor::letter a,b,c,i,m;

    libtensor::letter i,j,k,l,m,n,a,b,c,d,e,f;
    iw_ovvv(i|c|b|a) = i_ovvv(i|c|b|a)
                       - contract(m, t1(m|c), i_oovv(m|i|a|b));

        // iw_ovvv(k|a|c|d) =
        //     // 0.5 * asymm(c,d, 
        //       i_ovvv(k|a|c|d)
        //     - contract(l, i_oovv(k|l|c|d), t1(l|a))
        //     // )
        //     ;               
    return iw_ovvv;

}

#endif //GMB_EOM_CCSD_IW2_OVVV_H
