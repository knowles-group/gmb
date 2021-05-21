#ifndef GMB_EOM_CCSD_IW_OOOV_H
#define GMB_EOM_CCSD_IW_OOOV_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief W_ooov intermediate
 * 
 */
container<4,double> eom_ccsd_iw_ooov(
    container<2, double> &t1,       ///< CCSD T1
    container<4, double> &t2,       ///< CCSD T2 intermediate
    container<4, double> &tau,      ///< tau intermediate
    container<2, double> &f_ov,    ///< fock matrix - ov block
    container<2, double> &if_ov,    ///< fock matrix - ov block
    container<4, double> &iw_oooo,  ///> anti-symmetrized integral <ij||kl>
    container<4, double> &i_ooov,   ///> anti-symmetrized integral <ij||ka>
    container<4, double> &i_oovv,   ///> anti-symmetrized integral <ij||ab>
    container<4, double> &i_ovov,   ///> anti-symmetrized integral <ia||jb>
    container<4, double> &i_ovvv    ///> anti-symmetrized integral <ia||jb>

) {

    container<4, double> iw_ooov(i_ooov.get_space());

    // libtensor::letter a,i,j,k,m,e,f;

    libtensor::letter i,j,k,l,m,n,a,b,c,d,e,f;
    iw_ooov(i|j|k|a) = i_ooov(i|j|k|a)
                     - contract(e, t2(i|j|a|e), if_ov(k|e))
                     - contract(m, t1(m|a), iw_oooo(i|j|k|m))
                     + 0.5 * contract(e|f, tau(i|j|e|f), i_ovvv(k|a|e|f))
                     - asymm(i, j,
                        contract(e, t1(i|e), 
                            i_ovov(k|a|j|e)
                            - contract(m|f, t2(j|m|a|f), i_oovv(k|m|e|f)))
                       + contract(m|e, t2(i|m|a|e), i_ooov(k|m|j|e)));

        // iw_ooov(i|j|k|a) =
        //   // 0.5 * asymm(i,j, 
        //       i_ooov(i|j|k|a)
        //     - contract(l, iw_oooo(i|j|k|l), t1(l|a))
        //     + 0.5 * contract(c|d, i_ovvv(k|a|c|d), tau(i|j|c|d))
        //     - asymm(i, j,
        //         contract(b,
        //           i_ovov(k|a|j|b) - contract(l|c, i_oovv(k|l|b|c), t2(j|l|a|c)),
        //           t1(i|b))
        //       + contract(l|c, i_ooov(k|l|j|c), t2(i|l|a|c)))
        //     + contract(c, contract(l|b, i_oovv(k|l|b|c), t1(l|b)), t2(i|j|a|c))
        //     - contract(c, t2(i|j|a|c), f_ov(k|c))
        //     // )
            ;
    return iw_ooov;

}

#endif //GMB_EOM_CCSD_IW_OOOV_H
