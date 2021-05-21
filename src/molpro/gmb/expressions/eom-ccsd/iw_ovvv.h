#ifndef GMB_EOM_CCSD_IW_OVVV_H
#define GMB_EOM_CCSD_IW_OVVV_H

#include <libtensor/libtensor.h>
#include "../../container.h"

/**
 * @brief W_ovvv intermediate
 * 
 */
container<4,double> eom_ccsd_iw_ovvv(
    container<2, double> &t1,       ///< CCSD T1
    container<4, double> &t2,       ///< CCSD T2 amplitude
    container<4, double> &tau,      ///< tau intermediate
    container<2, double> &f_ov,    ///< fock matrix - ov block
    container<2, double> &if_ov,    ///< fock matrix - ov block
    container<4, double> &i_ooov,   ///> anti-symmetrized integral <ij||ka>
    container<4, double> &i_oovv,   ///> anti-symmetrized integral <ij||ab>
    container<4, double> &i_ovov,   ///> anti-symmetrized integral <ia||kb>
    container<4, double> &i_ovvv,   ///> anti-symmetrized integral <ia||jb>
    container<4, double> &iw_vvvv    ///> anti-symmetrized integral <ab||cd>

) {

    container<4, double> iw_ovvv(i_ovvv.get_space());

    // libtensor::letter a,b,c,i,m,n,e;

    libtensor::letter i,j,k,l,m,n,a,b,c,d,e,f;
    iw_ovvv(i|a|c|b) = i_ovvv(i|a|c|b)
                     - contract(m, t2(i|m|c|b), if_ov(m|a))
                     + contract(e, t1(i|e), iw_vvvv(b|c|a|e))
                     + asymm(b, c,
                        - contract(m, t1(m|c), i_ovov(m|b|i|a))
                        + contract(m|e, 
                            t2(i|m|c|e), i_ovvv(m|b|e|a))
                        - contract(m, t1(m|b), 
                            contract(e|n, t2(n|i|e|c), i_oovv(m|n|a|e))))
                     + 0.5 * contract(m|n, tau(m|n|c|b), i_ooov(m|n|i|a));
        // iw_ovvv(i|a|b|c) =
        //    0.5 * asymm(b,c, 
        //       i_ovvv(i|a|b|c)
        //     - contract(d, iw_vvvv(b|c|a|d), t1(i|d))
        //     + 0.5 * contract(j|k, i_ooov(j|k|i|a), t2(j|k|b|c))
        //     + asymm(b, c,
        //         contract(j,
        //             i_ovov(j|b|i|a)
        //           - 0.5 * contract(k, i_ooov(j|k|i|a), t1(k|b))
        //           - contract(k|d, i_oovv(j|k|a|d), t2(i|k|b|d)),
        //             t1(j|c))
        //       - contract(k|d, i_ovvv(k|c|a|d), t2(i|k|b|d)))
        //     - contract(k,
        //         contract(j|d, i_oovv(k|j|a|d), t1(j|d)),
        //         t2(i|k|b|c))
        //     - contract(k, t2(i|k|b|c), f_ov(k|a))
        //     )
        //     ;
    return iw_ovvv;

}

#endif //GMB_EOM_CCSD_IW_OVVV_H
