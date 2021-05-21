#ifndef GMB_EOM_CCSD_IW2_OOOV_H
#define GMB_EOM_CCSD_IW2_OOOV_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief W_ooov intermediate
 * 
 */
container<4,double> eom_ccsd_iw2_ooov(
    container<2, double> &t1,       ///< CCSD T1
    container<4, double> &i_ooov,   ///> anti-symmetrized integral <ij||ka>
    container<4, double> &i_oovv    ///> anti-symmetrized integral <ij||ab>

) {

    container<4, double> iw_ooov(i_ooov.get_space());

    libtensor::letter i,j,k,l,m,n,a,b,c,d,e,f;
    // libtensor::letter a,i,j,k,e;

    iw_ooov(j|k|i|a) = i_ooov(j|k|i|a)
                     + contract(e, t1(i|e), i_oovv(j|k|e|a));
        // iw_ooov(k|l|i|c) =
        // //   0.5 * asymm(k,l, 
        //       i_ooov(k|l|i|c)
        //     - contract(d, i_oovv(k|l|c|d), t1(i|d))
        //     // )
        //     ;
    return iw_ooov;

}

#endif //GMB_EOM_CCSD_IW2_OOOV_H
