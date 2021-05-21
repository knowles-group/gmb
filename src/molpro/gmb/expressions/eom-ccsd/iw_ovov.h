#ifndef GMB_CCSD_EOM_IW_OVOV_H
#define GMB_CCSD_EOM_IW_OVOV_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief W_ovov intermediate
 * 
 */
container<4,double> eom_ccsd_iw_ovov(
    container<2, double> &t1,       ///< CCSD T1 amplitude
    container<4, double> &t2,       ///< tau intermediate
    container<4, double> &i_ooov,   ///> anti-symmetrized integral <ij||ka>
    container<4, double> &i_oovv,   ///> anti-symmetrized integral <ij||ab>
    container<4, double> &i_ovov,   ///> anti-symmetrized integral <ia||jb>
    container<4, double> &i_ovvv    ///> anti-symmetrized integral <ia||bc>

) {

    container<4, double> iw_ovov(i_ovov.get_space());

    // libtensor::letter a,b,i,j,m,e;
    
    libtensor::letter i,j,k,l,m,n,a,b,c,d,e,f;
    iw_ovov(i|a|j|b) = i_ovov(i|a|j|b)
                     - contract(e, t1(i|e), i_ovvv(j|b|a|e))
                     - contract(m|e, t2(i|m|b|e), i_oovv(j|m|a|e))
                     + contract(m, 
                        contract(e, t1(i|e), i_oovv(j|m|a|e))
                         - i_ooov(j|m|i|a), t1(m|b));
        // iw_ovov(i|a|j|b) =
        //       i_ovov(i|a|j|b)
        //     - contract(c, i_ovvv(j|b|a|c), t1(i|c))
        //     - contract(k|c, i_oovv(j|k|a|c), t2(i|k|b|c))
        //     + contract(k,
        //         contract(c, i_oovv(j|k|a|c), t1(i|c)) - i_ooov(j|k|i|a),
        //         t1(k|b));
    return iw_ovov;

}

#endif //GMB_CCSD_EOM_IW_OVOV_H
